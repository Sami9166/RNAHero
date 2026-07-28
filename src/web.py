"""Local browser UI for RNAHero."""

from __future__ import annotations

import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from run import clear_output, run_disease


PAGE = """<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>RNAHero</title><style>
:root{color:#182230;background:#f7f8fb;font:16px/1.5 Inter,system-ui,sans-serif}*{box-sizing:border-box}body{margin:0}.wrap{width:min(920px,calc(100% - 32px));margin:0 auto;padding:72px 0}.brand{font-size:20px;font-weight:800;letter-spacing:-.05em}.brand span{color:#5b5ce2}h1{font-size:clamp(32px,6vw,52px);line-height:1.08;letter-spacing:-.06em;margin:24px 0 12px}.intro{max-width:620px;color:#657086;margin:0 0 28px}.ask{display:flex;gap:10px}.ask input{flex:1;min-width:0;border:1px solid #dce1eb;border-radius:14px;padding:16px 18px;font:inherit;background:#fff;box-shadow:0 8px 22px #1e2a4a0a}.ask button,.options button{border:0;border-radius:14px;padding:0 20px;font:inherit;font-weight:700;cursor:pointer}.ask button{color:#fff;background:#292c5f}.ask button:disabled{opacity:.6;cursor:wait}.options{display:flex;gap:10px;align-items:center;margin:14px 2px 34px;color:#657086;font-size:14px}.options input{width:120px;border:1px solid #dce1eb;border-radius:8px;padding:7px;background:#fff}.options button{padding:7px 10px;color:#657086;background:transparent}.panel{background:#fff;border:1px solid #e4e8f0;border-radius:20px;padding:26px;box-shadow:0 18px 50px #1e2a4a0a}.hidden{display:none}.status{display:inline-flex;align-items:center;background:#eef0ff;color:#45469d;border-radius:999px;padding:7px 11px;font-size:13px;font-weight:700}.status.error{background:#fff0f1;color:#a33c48}.timeline{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:24px 0 14px}.step{border-radius:12px;padding:12px;background:#f4f5f8;color:#8790a2;font-size:13px}.step.active{background:#e9e9ff;color:#41429c}.step.done{background:#edf8f2;color:#28704b}.step b{display:block;font-size:14px;margin-bottom:3px}.fact{min-height:24px;color:#657086;font-size:14px}.results{display:grid;grid-template-columns:repeat(2,1fr);gap:12px;margin-top:24px}.block{border:1px solid #e6e9f0;border-radius:14px;padding:16px}.block h3{font-size:14px;margin:0 0 10px}.metric{font-size:28px;font-weight:800;letter-spacing:-.05em}.muted{font-size:13px;color:#758095;margin:4px 0 0}.list{display:grid;gap:8px;margin-top:16px}.card{border:1px solid #e8ebf1;border-radius:12px;padding:12px}.tag{font-size:12px;color:#5b5ce2;font-weight:800;margin-bottom:4px}.card strong{display:block;font-size:14px}.card p{margin:4px 0 0}.score{float:right;color:#343675;font-weight:800}@media(max-width:620px){.wrap{padding:38px 0}.ask{display:grid}.ask button{height:52px}.timeline,.results{grid-template-columns:1fr 1fr}.options{flex-wrap:wrap}}
</style><body><main class="wrap"><div class="brand">RNA<span>Hero</span></div><h1>질환명을 입력하면<br>근거 기반 마커를 찾습니다.</h1><p class="intro">GEO raw count를 수집하고, sample metadata를 검토한 뒤 edgeR 내부·외부 검증까지 진행합니다.</p><form id="form" class="ask"><input id="disease" aria-label="질환명" placeholder="예: lung adenocarcinoma" required><button>분석 시작</button></form><div class="options">저장 폴더 <input id="output" value="output" aria-label="저장 폴더"><button id="clear" type="button">기록 지우기</button></div><section id="result" class="panel hidden"><div id="status" class="status"></div><div id="timeline" class="timeline"></div><p id="fact" class="fact"></p><div id="results" class="results"></div><div id="lists" class="list"></div></section></main><script>
const form=document.querySelector('#form'),button=form.querySelector('button'),clear=document.querySelector('#clear'),result=document.querySelector('#result'),status=document.querySelector('#status'),timeline=document.querySelector('#timeline'),fact=document.querySelector('#fact'),results=document.querySelector('#results'),lists=document.querySelector('#lists'),output=document.querySelector('#output');
const steps=[['searching','GEO 검색','raw-count 후보를 찾습니다'],['selecting_cohorts','표본 검토','case/control metadata를 확인합니다'],['analysing','edgeR 검증','세 코호트로 내부 검증합니다'],['complete','결과 정리','마커와 산출물을 저장합니다']];
const facts=['GEO 제목만으로 cohort를 확정하지 않습니다. sample metadata를 함께 검토합니다.','AUC는 두 내부 검증 코호트 중 더 낮은 값으로 보수적으로 평가합니다.','raw count는 보존하고, edgeR 결과와 검증 점수는 별도 파일로 저장합니다.'];let poller,factTimer,factIndex=0,activeOutput='output';
function stageIndex(stage){const index=steps.findIndex(step=>step[0]===stage);return index<0?0:index}function text(node,value){node.textContent=value;return node}function el(tag,cls){const node=document.createElement(tag);if(cls)node.className=cls;return node}
function drawTimeline(stage){const active=stageIndex(stage);timeline.innerHTML='';steps.forEach((step,index)=>{const state=index<active?'done':index===active&&stage!=='complete'?'active':stage==='complete'?'done':'';const card=el('div','step '+state);card.append(text(el('b'),step[1]),text(el('span'),step[2]));timeline.append(card)})}
function drawList(title,items,renderer){if(!items.length)return;const block=el('section','block');block.append(text(el('h3'),title));items.forEach(item=>block.append(renderer(item)));lists.append(block)}
function studyCard(study){const card=el('article','card');card.append(text(el('div','tag'),study.gse_id),text(el('strong'),study.title||study.gse_id));const detail=study.role?((study.role==='external'?'외부 검증':'개발')+' · case '+study.case_count+' / control '+study.control_count):'검색 후보';card.append(text(el('p','muted'),detail));return card}
function biomarkerCard(gene){const card=el('article','card');card.append(text(el('span','score'),'AUC '+Number(gene.min_auc).toFixed(3)),text(el('strong'),gene.gene_id),text(el('p','muted'),'민감도 '+Number(gene.min_sensitivity).toFixed(3)+' · 특이도 '+Number(gene.min_specificity).toFixed(3)));return card}
function criticCard(critic){const card=el('article','card');card.append(text(el('div','tag'),'ADK · '+(critic.verdict||critic.status||'unavailable')),text(el('strong'),'검증 점검'),text(el('p','muted'),critic.summary||critic.reason||'critic report unavailable'));return card}
function render(data){result.classList.remove('hidden');const running=['searching','selecting_cohorts','analysing'].includes(data.stage);status.className='status'+(data.stage.includes('failed')?' error':'');text(status,'상태: '+data.stage+(data.reason?' — '+data.reason:''));drawTimeline(data.stage);fact.hidden=!running;if(!running)fact.textContent='';lists.innerHTML='';results.innerHTML='';const analysis=data.analysis||{};if(data.stage==='complete'){[['최종 마커',analysis.biomarker_count||0],['점수화 후보',analysis.candidate_score_count||0]].forEach(pair=>{const block=el('section','block');block.append(text(el('h3'),pair[0]),text(el('div','metric'),String(pair[1])));results.append(block)});drawList('ADK critic 점검',data.critic?[data.critic]:[],criticCard);drawList('실제 분석 코호트',data.cohorts||[],studyCard);drawList('상위 바이오마커',analysis.top_biomarkers||[],biomarkerCard)}else{drawList('검색 후보',data.studies||[],studyCard)}}
async function latest(){try{const response=await fetch('/api/latest?output='+encodeURIComponent(activeOutput));if(response.ok)render(await response.json())}catch(_){}}
function startPolling(){stopPolling();latest();poller=setInterval(latest,1500);factTimer=setInterval(()=>{factIndex=(factIndex+1)%facts.length;fact.textContent=facts[factIndex]},4200);fact.textContent=facts[0]}function stopPolling(){clearInterval(poller);clearInterval(factTimer)}
form.addEventListener('submit',async event=>{event.preventDefault();activeOutput=output.value.trim()||'output';button.disabled=true;button.textContent='분석 중…';startPolling();try{const response=await fetch('/api/run',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({disease:document.querySelector('#disease').value,output:activeOutput})});render(await response.json())}catch(error){render({stage:'failed',reason:error.message,studies:[]})}finally{stopPolling();button.disabled=false;button.textContent='분석 시작'}});
clear.addEventListener('click',async()=>{if(!confirm('현재 저장 폴더의 검색·분석 기록을 지울까요?'))return;activeOutput=output.value.trim()||'output';const response=await fetch('/api/clear',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({output:activeOutput})});const data=await response.json();render(data.stage==='cleared'?{stage:'ready',studies:[]}:data)});activeOutput=output.value;latest();
</script></body></html>"""


class Handler(BaseHTTPRequestHandler):
    def _json(self, value: object, status: int = 200) -> None:
        data = json.dumps(value, ensure_ascii=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/":
            data = PAGE.encode()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        elif parsed.path == "/api/latest":
            output = parse_qs(parsed.query).get("output", ["output"])[0]
            path = Path(output) / "provenance" / "run_manifest.json"
            if path.is_file():
                self._json(json.loads(path.read_text(encoding="utf-8")))
            else:
                self._json({"stage": "ready", "studies": []})
        else:
            self._json({"error": "not found"}, 404)

    def do_POST(self) -> None:
        try:
            size = int(self.headers.get("Content-Length", "0"))
            request = json.loads(self.rfile.read(size))
            if self.path == "/api/clear":
                self._json(clear_output(Path(str(request.get("output", "output")))))
                return
            if self.path != "/api/run":
                self._json({"error": "not found"}, 404)
                return
            disease = str(request["disease"]).strip()
            if not disease:
                raise ValueError("질환명을 입력하세요")
            self._json(run_disease(disease, Path(str(request.get("output", "output")))))
        except Exception as error:
            self._json({"stage": "failed", "reason": str(error), "studies": []}, 400)

    def log_message(self, _format: str, *_args: object) -> None:
        return


def serve(port: int = 8000) -> None:
    server = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    print(f"RNAHero is running at http://127.0.0.1:{port}")
    server.serve_forever()
