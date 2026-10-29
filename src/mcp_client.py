"""Synchronous bridge from the RNAHero pipeline to the local MCP server."""

from __future__ import annotations

import asyncio
import atexit
import json
import sys
import threading
from concurrent.futures import Future
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parent.parent
_CLOSE_TIMEOUT_SECONDS = 5
_START_TIMEOUT_SECONDS = 60

_lock = threading.Lock()
_loop: asyncio.AbstractEventLoop | None = None
_session: Any = None
_stop: asyncio.Event | None = None


async def _serve(ready: Future, stopped: Future) -> None:
    """Own the stdio session for its whole life so anyio scopes stay in one task."""
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    global _stop
    try:
        parameters = StdioServerParameters(command=sys.executable, args=["-m", "ncbi_mcp"], cwd=str(PROJECT_ROOT))
        async with stdio_client(parameters) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                _stop = asyncio.Event()
                ready.set_result(session)
                await _stop.wait()
    except Exception as error:
        if not ready.done():
            ready.set_exception(error)
    finally:
        stopped.set_result(None)


def _ensure() -> tuple[asyncio.AbstractEventLoop, Any]:
    global _loop, _session
    with _lock:
        if _session is None:
            loop = asyncio.new_event_loop()
            threading.Thread(target=loop.run_forever, name="rnahero-mcp", daemon=True).start()
            ready: Future = Future()
            stopped: Future = Future()
            asyncio.run_coroutine_threadsafe(_serve(ready, stopped), loop)
            try:
                _session = ready.result(timeout=_START_TIMEOUT_SECONDS)
            except Exception:
                loop.call_soon_threadsafe(loop.stop)
                raise
            _loop, _stopped = loop, stopped
            atexit.register(lambda: close(_stopped))
        return _loop, _session


def call_tool(name: str, **arguments: Any) -> Any:
    """Call one tool on the local MCP server and return its decoded JSON payload."""
    loop, session = _ensure()
    result = asyncio.run_coroutine_threadsafe(session.call_tool(name, arguments), loop).result()
    text = "".join(item.text for item in result.content if getattr(item, "text", None))
    if result.isError:
        raise RuntimeError(f"MCP tool {name} failed: {text or result}")
    if text:
        return json.loads(text)
    # An empty list or dict arrives only as structured output, which FastMCP wraps for non-objects.
    structured = getattr(result, "structuredContent", None) or {}
    return structured.get("result", structured)


def close(stopped: Future | None = None) -> None:
    """Stop the server subprocess; a hung shutdown must not block the caller."""
    global _loop, _session, _stop
    with _lock:
        loop, stop = _loop, _stop
        _loop, _session, _stop = None, None, None
    if loop is None:
        return
    if stop is not None:
        loop.call_soon_threadsafe(stop.set)
    if stopped is not None:
        try:
            stopped.result(timeout=_CLOSE_TIMEOUT_SECONDS)
        except Exception:
            pass
    loop.call_soon_threadsafe(loop.stop)
