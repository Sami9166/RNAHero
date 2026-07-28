# RNAHero

GEO의 bulk RNA-seq raw count로 연구용 바이오마커 후보를 찾는 로컬 파이프라인입니다. 질환명(한글·영문)을 입력하면 GEO 검색부터 cohort 판정, edgeR, 내부·외부 검증, 비판적 검토와 Markdown 보고서까지 실행합니다.

> 결과는 연구 후보이며 임상 진단 또는 치료 의사결정에 사용할 수 없습니다.

## 설치

```powershell
cd C:\Users\Min\Documents\RNAHero
python -m venv .venv
.\.venv\Scripts\python -m pip install -e ".[agent]"
```

`.env`에 Gemini 키를 설정합니다.

```text
GEMINI_API_KEY=your_key
GEMINI_MODEL=gemini-3.5-flash-lite
```

선택적으로 NCBI E-utilities 검색 한도를 늘릴 수 있습니다.

```text
NCBI_API_KEY=your_ncbi_key
NCBI_EMAIL=you@example.com
```

## 빠른 실행

```powershell
.\.venv\Scripts\rnahero.exe run "lung adenocarcinoma" --output output
.\.venv\Scripts\rnahero.exe web
```

브라우저 UI는 `http://127.0.0.1:8000`에서 열립니다. 중단된 분석은 이미 저장된 GEO 파일과 완료된 edgeR 결과를 재사용합니다.

```powershell
.\.venv\Scripts\rnahero.exe resume --output output
```

완료된 분석의 보고서와 그림만 다시 만들려면 다음을 실행합니다.

```powershell
.\.venv\Scripts\rnahero.exe summarize --output output
```

## 파이프라인

```text
질환명
  ↓
NCBI GEO MCP: bulk RNA-seq + raw-count 후보 최소 20개 검색
  ↓
Gemini ADK: 각 GSM metadata에서 환자 case/control 근거 확인
  └─ 세포주·약물/유전자 처리·비임상 검체·근거 충돌은 제외
  ↓
3개 development cohort + 1개 held-out external cohort 선택
  ↓
raw count 열 ↔ GSM ↔ case/control crosswalk 생성 및 edgeR 입력 정리
  ↓
development 3개 edgeR 병렬 실행
  ↓
각 discovery cohort의 FDR 후보를 나머지 development 2개에서 내부 검증
  └─ AUC ≥ 0.80, 민감도 ≥ 0.70, 특이도 ≥ 0.70을 모두 충족한 후보만 lock
  ↓
lock된 후보만 external cohort에서 평가
  ↓
ADK critic: cohort 근거·내부/외부 검증·재현성 한계 비판적 점검
  ↓
GO enrichment: primary development cohort의 FDR 통과 DEG 기능 분석과 dot plot 생성
  ↓
ADK summarizer: 상위 5개 logFC 후보, GO 결과, critic 비판점, PCA·heatmap 경로를 Markdown 보고서로 작성
```

external cohort는 후보를 고르거나 기준을 바꾸는 데 사용하지 않습니다. 후보 순위는 development cohort에서 내부 검증을 통과하고 FDR 기준을 만족한 유전자들의 절대 logFC를 기준으로 정합니다.

## Agent 역할

| 구성 요소 | 역할 | 변경하지 않는 것 |
| --- | --- | --- |
| NCBI MCP | GEO 검색, GSM metadata, raw-count 다운로드 | 바이오마커 결론 |
| Gemini orchestrator | sample field 근거를 바탕으로 cohort 적합성 판정 | edgeR 결과와 cut-off |
| edgeR | 차등발현과 log-CPM 산출 | cohort 역할 |
| ADK critic | 전체 파이프라인과 검증 결과의 한계 점검 | 후보·방향·FDR·기준 |
| ADK summarizer | 고정된 상위 5개 후보와 그림을 사람이 읽는 보고서로 정리 | 후보 재순위화 |

## 산출물

```text
output/
  search_results.json                 # GEO 후보 목록
  cohorts.json                        # 실제 분석에 채택된 4개 cohort
  metadata/GSE....json                # 채택 cohort의 GEO metadata
  raw/GSE.../                         # 내려받은 원본 raw count
  prepared/GSE....csv                 # edgeR 입력 행렬
  samples/GSE....csv                  # case/control sample sheet
  analysis/
    analysis_config.json
    GSE.../crosswalk.json             # count 열 ↔ GSM ↔ group
    GSE.../edger_results.csv
    GSE.../logcpm.csv
    biomarkers.json / biomarkers.csv
    internal_validation.json / .csv
    external_validation.json / .csv
    validation_report.json
    critic_report.json
    summary_report.md                 # ADK summarizer의 최종 보고서
    summary/figures/GSE.../
      top5_pca.png
      top5_heatmap.png
    summary/go/GSE.../
      go_enrichment.csv                # FDR 통과 DEG의 GO enrichment 결과
      go_dotplot.png
```

PCA와 heatmap은 각각 development cohort와 external cohort에 생성됩니다. PCA는 샘플 수준의 군 분리를, heatmap은 상위 5개 유전자의 발현 패턴을 보여 줍니다.
GO enrichment은 primary development cohort의 FDR 통과 DEG를 up/down으로 나누고 g:Profiler의 인간 유전체 기본 배경으로 수행합니다. 따라서 기능 해석은 cohort-specific이며, 외부 cohort로 후보를 선택하지 않습니다.

## 수동 edgeR 실행

count 파일은 첫 열이 `gene_id`이고 나머지가 정수 raw count여야 합니다. sample sheet는 `sample_id,group` 형식이며 group은 `case` 또는 `control`입니다. 다음처럼 development 3개와 external 1개를 둔 config를 작성할 수 있습니다.

```json
{
  "fdr": 0.05,
  "output_dir": "outputs",
  "development": [
    {"name": "dev1", "counts": "data/dev1_counts.csv", "samples": "data/dev1_samples.csv"},
    {"name": "dev2", "counts": "data/dev2_counts.csv", "samples": "data/dev2_samples.csv"},
    {"name": "dev3", "counts": "data/dev3_counts.csv", "samples": "data/dev3_samples.csv"}
  ],
  "external": {"name": "external", "counts": "data/external_counts.csv", "samples": "data/external_samples.csv"}
}
```

```powershell
.\.venv\Scripts\rnahero.exe pipeline config.json
```

## 한계

- GEO metadata만으로 적합한 환자-대조군 cohort가 충분하지 않으면 분석은 `needs_cohort_review`에서 멈춥니다.
- Gene symbol과 Ensembl ID는 로컬 캐시와 Ensembl fallback으로 표준화하며, 매핑할 수 없는 유전자는 검증에서 제외됩니다.
- 코호트 간 log-CPM cut-off 비교는 MVP 방식입니다. 임상적 주장을 하려면 배치 보정, 사전등록된 모델·cut-off, 더 큰 독립 검증이 필요합니다.

## 테스트

```powershell
$env:PYTHONPATH = "src"
.\.venv\Scripts\python -m unittest discover -s tests -v
```
