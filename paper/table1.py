"""Recompute manuscript Table 1 from the released evaluation output (no API calls).

    python paper/table1.py
"""

import json
from pathlib import Path

from pathwayseeker.evaluation import summarize

HERE = Path(__file__).parent


def category(r):
    qid = r["query_id"]
    if qid.startswith("t1_connected"):
        return "Connected"
    if qid.startswith("t1_unconnected"):
        return "Unconnected"
    return "Phenylpropanoid"


def main():
    data = json.loads((HERE / "results" / "table1_results.json").read_text())
    rows = data["by_mode"]["hypothesis"]["results"]
    for r in rows:
        r["eer"] = r.get("support_ratio", 0.0)
        r["overall"] = r.get("quality_overall")
    table = summarize(rows, category_of=category)
    cols = ["n", "eer_pct", "scientific_reasoning", "specificity", "evidence_transparency", "clarity", "overall"]
    print(f"{'Category':<16}" + "".join(f"{c:>22}" for c in cols))
    for cat in ("Connected", "Unconnected", "Phenylpropanoid", "All queries"):
        print(f"{cat:<16}" + "".join(f"{str(table[cat][c]):>22}" for c in cols))
    return table


if __name__ == "__main__":
    main()
