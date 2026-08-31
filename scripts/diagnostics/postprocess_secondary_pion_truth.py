#!/usr/bin/env python3
"""Make plots and accounting from the compact secondary-pion truth ROOT file."""
import argparse, collections, json, time
from pathlib import Path
import ROOT
ROOT.gROOT.SetBatch(True)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("--output-dir",required=True); a=ap.parse_args()
    start=time.time(); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True); f=ROOT.TFile.Open(a.input); t=f.Get("EventSummary"); d=f.Get("DirectDaughterPi")
    rows=[]
    for x in t:
        rows.append({n:int(getattr(x,n)) for n in ("PrimaryReachedS2","HasDirectDaughterPi","HasDirectDaughterReachedHGC","HasDirectDaughterReachedAGC","HasDirectDaughterReachedS2","HasDirectDaughterReachedCal","HasDirectDaughterExitedCal","DirectDaughterPiInPIDRegion","DescendantPiInPIDRegion","AnyPiInPIDRegion","TriggerPass","OldPionPIDPass")}|{"AGCNPE":int(x.AGCNPE),"HGCNPE":int(x.HGCNPE)})
    n=len(rows); missed=[r for r in rows if not r["PrimaryReachedS2"]]; mt=[r for r in missed if r["TriggerPass"]]; mtp=[r for r in mt if r["OldPionPIDPass"]]
    def count(rs,key): return sum(r[key] for r in rs)
    # Exclusive categories: direct takes precedence, then later descendant, then no pion.
    def cat(r):
        if r["DirectDaughterPiInPIDRegion"]: return "direct daughter"
        if r["DescendantPiInPIDRegion"]: return "later descendant"
        if not r["AnyPiInPIDRegion"]: return "no pi+ in PID"
        return "other/unrelated pi+"
    table={c:{"pass":0,"fail":0} for c in ("direct daughter","later descendant","no pi+ in PID","other/unrelated pi+")}
    for r in mt: table[cat(r)]["pass" if r["OldPionPIDPass"] else "fail"]+=1
    summary={"generated":n,"primary_reached_s2":count(rows,"PrimaryReachedS2"),"primary_missed_s2":len(missed),"primary_missed":{"direct_produced":count(missed,"HasDirectDaughterPi"),"direct_reached_hgc":count(missed,"HasDirectDaughterReachedHGC"),"direct_reached_agc":count(missed,"HasDirectDaughterReachedAGC"),"direct_reached_s2":count(missed,"HasDirectDaughterReachedS2"),"direct_reached_cal":count(missed,"HasDirectDaughterReachedCal"),"direct_exited_cal":count(missed,"HasDirectDaughterExitedCal")},"primary_missed_trigger_pass":len(mt),"missed_trigger_pass":{"direct_in_pid":count(mt,"DirectDaughterPiInPIDRegion"),"descendant_in_pid":count(mt,"DescendantPiInPIDRegion"),"no_pi_in_pid":sum(not r["AnyPiInPIDRegion"] for r in mt)},"primary_missed_trigger_pion_pid_pass":len(mtp),"two_dimensional_accounting":table,"timing_postprocess_s":time.time()-start}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    cats={"primary_survives":"PrimaryReachedS2","direct_daughter_recovery":"!PrimaryReachedS2&&DirectDaughterPiInPIDRegion","later_descendant":"!PrimaryReachedS2&&!DirectDaughterPiInPIDRegion&&DescendantPiInPIDRegion","no_pion_in_PID":"!PrimaryReachedS2&&!AnyPiInPIDRegion"}
    colors=[ROOT.kBlue+1,ROOT.kGreen+2,ROOT.kOrange+7,ROOT.kRed+1]
    for var,bins,lo,hi in (("AGCNPE",120,0,600),("HGCNPE",120,0,3000)):
        c=ROOT.TCanvas("c","c",1000,750); leg=ROOT.TLegend(.58,.62,.88,.88); hs=[]
        for i,(label,cut) in enumerate(cats.items()):
            h=ROOT.TH1D(f"h_{var}_{i}",f"{var};{var};normalized events",bins,lo,hi);t.Draw(f"{var}>>{h.GetName()}",cut,"goff");h.SetDirectory(0);h.SetLineColor(colors[i]);h.SetLineWidth(3); h.Scale(1/h.Integral() if h.Integral() else 1);hs.append(h);leg.AddEntry(h,label.replace('_',' '),"l")
        ymax=max(h.GetMaximum() for h in hs)*1.25;hs[0].SetMaximum(ymax);hs[0].Draw("hist");[h.Draw("hist same") for h in hs[1:]];line=ROOT.TLine(20 if var=="AGCNPE" else 50,0,20 if var=="AGCNPE" else 50,ymax);line.SetLineStyle(2);line.Draw();leg.Draw();c.SaveAs(str(out/f"{var}_truth_categories.png"))
    for label,cut in cats.items():
        c=ROOT.TCanvas("c2","c2",900,750);h=ROOT.TH2D("h2",f"{label.replace('_',' ')};AGCNPE;HGCNPE",100,0,600,100,0,3000);t.Draw("HGCNPE:AGCNPE>>h2",cut,"colz");ROOT.TLine(20,0,20,3000).Draw();ROOT.TLine(0,50,600,50).Draw();c.SaveAs(str(out/f"AGC_vs_HGC_{label}.png"));h.Delete()
    c=ROOT.TCanvas("c3","c3",1000,750);h=ROOT.TH1D("stage","Direct daughter furthest stage;stage;tracks",7,-.5,6.5);d.Draw("FurthestStage>>stage","","goff");[h.GetXaxis().SetBinLabel(i+1,x) for i,x in enumerate(("BeforeS1","S1","HGC","AGC","S2","Cal","ExitedCal"))];h.Draw("hist");c.SaveAs(str(out/"direct_daughter_furthest_stage.png"))
    pct=lambda x,den:100*x/den if den else 0
    md=["# Secondary-pion visibility study","","> AGCNPE and HGCNPE are event-level quantities. Truth-pion presence is correlated with these responses; the file cannot prove that a particular track caused the NPE.","","`InPIDRegion` is defined as `ReachedHGC OR ReachedAGC`, following the actual S1 → HGC → AGC → S2 → Cal geometry.","",f"Generated events: **{n:,}**",f"Primary pi+ reached S2: **{summary['primary_reached_s2']:,} ({pct(summary['primary_reached_s2'],n):.3f}%)**",f"Primary pi+ missed S2: **{len(missed):,} ({pct(len(missed),n):.3f}%)**","","## Among primary-missed-S2 events",""]
    for k,v in summary["primary_missed"].items():md.append(f"- {k.replace('_',' ')}: {v:,} ({pct(v,len(missed)):.3f}%)")
    md += ["","## Trigger and old PID cross-check","",f"Primary missed S2 + trigger pass: **{len(mt):,}**",f"Primary missed S2 + trigger pass + old pion PID: **{len(mtp):,}**","","| Exclusive truth category | Old pion PID pass | Old pion PID fail |","|---|---:|---:|"]
    for k,v in table.items():md.append(f"| {k} | {v['pass']:,} | {v['fail']:,} |")
    md += ["","StartP and EndP are available for direct daughters. No momentum-at-HGC/AGC-crossing branch exists, so these are not presented as detector-crossing momenta."]
    (out/"SECONDARY_PION_VISIBILITY_STUDY.md").write_text("\n".join(md)+"\n")
    print(json.dumps(summary,indent=2));print(f"TIMING plotting_report_s {time.time()-start:.3f}")
if __name__=="__main__": main()
