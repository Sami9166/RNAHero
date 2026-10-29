import unittest
from concurrent.futures import ThreadPoolExecutor

import mcp_client


PIPELINE_TOOLS = {
    "search_bulk_rnaseq",
    "get_geo_sample_metadata",
    "match_geo_sample_columns",
    "download_raw_count_matrix",
    "create_pca_heatmap",
    "create_go_enrichment",
}
TOOL_TIMEOUT_SECONDS = 90


class McpClientTests(unittest.TestCase):
    """Exercise the real stdio server; no network tool is called."""

    @classmethod
    def tearDownClass(cls) -> None:
        mcp_client.close()

    def test_server_exposes_every_tool_the_pipeline_calls(self) -> None:
        import asyncio

        loop, session = mcp_client._ensure()
        tools = asyncio.run_coroutine_threadsafe(session.list_tools(), loop).result(TOOL_TIMEOUT_SECONDS)
        self.assertLessEqual(PIPELINE_TOOLS, {tool.name for tool in tools.tools})

    def test_tool_error_is_raised_without_blocking_the_server(self) -> None:
        # A lazy heavy import inside a tool used to deadlock the stdio loop, so this must return fast.
        with ThreadPoolExecutor(max_workers=1) as executor:
            call = executor.submit(
                mcp_client.call_tool, "create_pca_heatmap", output_dir="/outside", cohort_id="GSE1", gene_ids=["A", "B"]
            )
            with self.assertRaises(RuntimeError) as raised:
                call.result(TOOL_TIMEOUT_SECONDS)
        self.assertIn("create_pca_heatmap", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
