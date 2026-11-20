(* extrai_viga.wls: le os .mx da viga (Dropbox/UnsaturatedSoil/mathematica) e grava viga_dados.json, que o
   gera_dados.py le. Uso: wolframscript -file extrai_viga.wls   (leva ~1 min; so precisa rodar de novo se os .mx mudarem) *)
dir = FileNameJoin[{$HomeDirectory, "Dropbox", "UnsaturatedSoil", "mathematica"}];
ld[f_] := Import[FileNameJoin[{dir, f}]];
num[x_] := N[x];
ai[f_] := Module[{d = ld[f]}, <|"arquivo" -> f, "M" -> d["M"], "Mform" -> d["Mform"], "flechas" -> num[d["flechas"]],
    "pesos" -> num[d["pesos"]], "tempo" -> num[d["tempo"]], "analises" -> d["analises"]|>];
lsx[f_] := Module[{d = ld[f]}, <|"arquivo" -> f, "c" -> num[d["c"]], "nLinha" -> d["nLinha"], "s" -> num[d["s"]],
    "beta" -> num[d["beta"]], "tempo" -> num[d["tempo"]]|>];
fo[f_] := Module[{d = ld[f]}, <|"arquivo" -> f, "M" -> d["M"], "beta" -> num[d["beta"]], "tempo" -> num[d["tempo"]],
    "analises" -> d["analises"]|>];
mc = ld["mc_M200_N50000.mx"];
out = <|
  "ai" -> Join[
     <|"f8" -> ai["ai_M200_f8_N1000_d0.00075_s2.mx"], "f20" -> ai["ai_M200_f20_N1000_d0.00075_s2.mx"],
       "f200_s2" -> ai["ai_M200_f200_N1000_d0.00075_s2.mx"]|>,
     Association[Table["f200_s" <> ToString[k] -> ai["ai_M200_f200_N1000_d0.00075_s" <> ToString[k] <> ".mx"],
       {k, 101, 110}]]],
  "ls" -> <|"f8" -> lsx["ls_M200_f8_N500_d0.00075_s2.mx"], "f20" -> lsx["ls_M200_f20_N500_d0.00075_s2.mx"],
     "f200" -> lsx["ls_M200_f200_N500_d0.00075_s2.mx"], "adj200" -> lsx["ls_M200_fadj200_N500_d0.00075_s2.mx"]|>,
  "form" -> <|"M8" -> fo["form_M8_d0.00075.mx"], "M20" -> fo["form_M20_d0.00075.mx"],
     "M200" -> fo["form_M200_d0.00075.mx"], "adj200" -> fo["form_adj_M200_d0.00075.mx"]|>,
  "mc_flechas" -> num[mc[[2]]]|>;
Export[FileNameJoin[{DirectoryName[$InputFileName], "viga_dados.json"}], out, "RawJSON"];
Print["ok"];
