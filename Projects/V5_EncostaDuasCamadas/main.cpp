// V5_EncostaDuasCamadas: o caso completo de Abbasov et al. (2024) (ver README.md e problema.h). Etapa 1: ANOS anos de
// chuva a CHUVA1_MM_ANO a partir do estado hidrostatico com o nivel d'agua em z = NIVEL (ou o estado lido de
// ESTADO_INICIAL); etapa 2: CHUVA2_MM_H por HORAS horas; em cada instante de SAIDAS_H (horas da etapa 2; 0 = o fim da
// etapa 1) o calculo de tensoes elastico e o LFS na zona de investigacao.
//
// Ambiente (padroes entre parenteses; o resto no README.md): CASO (tabela1 | tabela4 | fig9), MALHA (malhas/encosta.msh),
// NU, PHI_DEG, C_SUP, C_INF (ou COES para as duas), KS_INF (cm/h), DT0_H, FATOR_DT, DT_MAX_D, DT2_H, SAIDAS_H
// (0,7.5,15), ESTADO_INICIAL, GRAVA_ESTADO, VTK (0), PONTOS_CSV (0), CAMPOS_CSV (1), VERBOSO (1), SOLVER (pardiso),
// NU_LISTA, C_LISTA, PHI_LISTA (sensibilidade das metricas do LFS em cada instante, em sensibilidade.csv).
//
// Saidas na pasta corrente: passos.csv (cada passo aceito, com o balanco de agua), metricas.csv (LFS medio, area com
// LFS < 1 e umidades de cada instante), perfis.csv (p_w, theta, S_e, tensao de succao e LFS em quatro verticais),
// campos_t<h>.csv (por elemento, para os mapas), nos.csv e elementos.csv (a malha), estado_etapa1.csv (p_w nos nos no
// fim da etapa 1), resumo.csv (uma linha) e resultado.txt (chave=valor, gravado por ultimo, de um .tmp renomeado).
// Codigo de saida 0 em sucesso.
#include "problema.h"

#include <cstdio>

int main() {
    using namespace v5;
    const auto inicio = Relogio::now();
    std::string situacao = "ok";
    try {
        const Config cf = Config::DoAmbiente();
        Verboso() = EnvI("VERBOSO", 1);
        TPZUnsatAnalysis::SetUseSkyline(EnvS("SOLVER", "pardiso") == "skyline");
        const std::string arqMalha = CaminhoMalha(cf.malha);
        Encosta e(cf, LeMalha(arqMalha));
        int nel = 0;
        for (TPZGeoEl *g : e.gmesh->ElementVec()) nel += g && g->Dimension() == 2;
        if (Verboso())
            std::cout << "V5 " << cf.caso << ": malha " << arqMalha << " (" << e.gmesh->NNodes() << " nos, " << nel
                      << " triangulos, " << e.an->NEq() << " equacoes, " << e.pontos.size() << " pontos), area da zona "
                      << e.areaZona << " m2\n  nu = " << cf.nu << ", phi' = " << cf.phiDeg << ", c' = " << cf.sup.c << " / "
                      << cf.inf.c << " kPa, K_s = " << cf.sup.ksCmH << " / " << cf.inf.ksCmH << " cm/h, VG a = "
                      << cf.sup.aVG << " / " << cf.inf.aVG << " 1/kPa, b = " << cf.sup.bVG << " / " << cf.inf.bVG
                      << std::endl;
        e.GravaMalha();
        std::ofstream fpas("passos.csv");
        fpas << std::setprecision(10);
        Encosta::CabecalhoPassos(fpas);

        // ---- etapa 1: 10 anos a 600 mm/ano (ou o estado lido)
        bool ok = true;
        const auto t1 = Relogio::now();
        int passos1 = 0;
        if (!cf.estadoIni.empty()) {
            e.LeEstado(cf.estadoIni);
            if (Verboso()) std::cout << "etapa 1: estado lido de " << cf.estadoIni << std::endl;
        } else {
            e.CondicaoInicial();
            ok = e.Etapa1(&fpas);
            passos1 = e.passos;
            if (Verboso())
                std::cout << "etapa 1: " << (ok ? "ok" : "NAO CONVERGIU") << ", " << e.passos << " passos, " << e.cortes
                          << " cortes, " << e.newtonTotal << " iteracoes, " << Segundos(t1) << " s; ultimo dt = "
                          << e.ultimo.dt / kDia << " d, max|dp_w| = " << e.ultimo.dpwMax << " kPa, face = "
                          << e.ultimo.faceNos << " nos" << std::endl;
            if (!ok) situacao = "etapa1_nao_convergiu";
        }
        // balanco do ultimo passo da etapa 1 (m3/s por metro): perto do permanente, chuva = - (dir + face)
        const InfoPasso fim1 = e.ultimo;
        e.GravaEstado("estado_etapa1.csv");
        if (!cf.gravaEstado.empty()) e.GravaEstado(cf.gravaEstado);

        // ---- etapa 2: 20 mm/h por 15 h; tensoes e LFS em cada instante
        std::ofstream fmet("metricas.csv"), fper("perfis.csv");
        fmet << "caso,t_h,LFS_medio,A_LFS_m2,area_zona_m2,LFS_min,LFS_medio_lim10,LFS_medio_lim3,A_LFS_sup_m2,A_LFS_inf_m2,"
                "A_LFS_maior3_m2,"
                "theta_sup_min,theta_sup_max,theta_inf_min,theta_inf_max,sigs_max_kPa,nivel_crista_m,frente_crista_m,"
                "frente_meio_m,theta_sup_crista,theta_inf_crista,sigs_sup_crista_kPa,sigs_inf_crista_kPa,theta_sup_pe,"
                "face_nos\n"
             << std::setprecision(8);
        // sensibilidade (opcional): NU_LISTA, C_LISTA (kPa, as duas camadas) e PHI_LISTA (graus) em cada instante
        const std::vector<double> nus = EnvLista("NU_LISTA", {}), cs = EnvLista("C_LISTA", {}),
                                  phis = EnvLista("PHI_LISTA", {});
        std::ofstream fsen;
        if (!nus.empty() || !cs.empty() || !phis.empty()) {
            fsen.open("sensibilidade.csv");
            fsen << "t_h,parametro,valor,LFS_medio,LFS_medio_lim3,A_LFS_m2\n" << std::setprecision(8);
        }
        std::vector<Instante> inst;
        bool primeiro = true, tensoesOk = true;
        const int passos0 = e.passos;
        const auto t2 = Relogio::now();
        auto saida = [&](double th) {
            if (!e.Tensoes()) { tensoesOk = false; return; }
            const Instante r = AvaliaInstante(e, th);
            inst.push_back(r);
            const Metricas &m = r.m;
            fmet << cf.caso << ',' << th << ',' << m.lfsMedio << ',' << m.aLFS << ',' << m.areaZona << ',' << m.lfsMin << ','
                 << m.lfsMedioLim << ',' << m.lfsMedioLim3 << ',' << m.aLFSsup << ',' << m.aLFSinf << ',' << m.aLFS3 << ',' << m.thetaSupMin << ','
                 << m.thetaSupMax << ',' << m.thetaInfMin << ',' << m.thetaInfMax << ',' << m.sigsMax << ','
                 << r.nivelCrista << ',' << r.frenteCrista << ',' << r.frenteMeio << ',' << r.thetaSupCrista << ','
                 << r.thetaInfCrista << ',' << r.sigsSupCrista << ',' << r.sigsInfCrista << ',' << r.thetaPe << ','
                 << e.face.size() << '\n';
            fmet.flush();
            GravaPerfis(e, th, fper, primeiro);
            fper.flush();
            primeiro = false;
            const std::string sufixo = "_t" + Num(th, 6);
            if (cf.csvCampos) e.GravaCampos("campos" + sufixo + ".csv");
            if (cf.csvPontos) e.GravaPontos("pontos" + sufixo + ".csv");
            if (cf.vtk) e.GravaVTK("encosta" + sufixo + ".vtk");
            if (fsen.is_open() && !Sensibilidade(e, th, nus, cs, phis, fsen)) tensoesOk = false;
            if (Verboso())
                std::cout << "t = " << th << " h: LFS medio = " << m.lfsMedio << ", A(LFS < 1) = " << m.aLFS
                          << " m2 (zona " << m.areaZona << " m2), LFS min = " << m.lfsMin << "; nivel na crista z = "
                          << r.nivelCrista << " m, frente na crista " << r.frenteCrista << " m; theta sup "
                          << m.thetaSupMin << "-" << m.thetaSupMax << ", inf " << m.thetaInfMin << "-" << m.thetaInfMax
                          << "; face " << e.face.size() << " nos" << std::endl;
        };
        if (ok) {
            ok = e.Etapa2(saida, &fpas);
            if (!ok) situacao = "etapa2_nao_convergiu";
            else if (!tensoesOk) situacao = "tensoes_nao_convergiram";
        }
        if (Verboso())
            std::cout << "etapa 2: " << e.passos - passos0 << " passos, " << Segundos(t2) << " s; total "
                      << Segundos(inicio) << " s" << std::endl;

        // ---- resumo.csv e resultado.txt
        auto valor = [&](double th, auto f) -> std::string {
            for (const auto &r : inst)
                if (std::fabs(r.m.t_h - th) < 1.e-9) return Num(f(r), 8);
            return "";
        };
        std::ofstream fres("resumo.csv");
        fres << "caso,malha,nos,triangulos,equacoes,nu,phi_deg,c_sup,c_inf,ks_sup_cm_h,ks_inf_cm_h,passos_etapa1,"
                "passos_etapa2,newton_total,cortes,chuva1_m3_s,dir1_m3_s,face1_m3_s";
        std::ostringstream lin;
        lin << std::setprecision(10) << cf.caso << ',' << std::filesystem::path(arqMalha).filename().string() << ','
            << e.gmesh->NNodes() << ',' << nel << ',' << e.an->NEq() << ',' << cf.nu << ',' << cf.phiDeg << ',' << cf.sup.c
            << ',' << cf.inf.c << ',' << cf.sup.ksCmH << ',' << cf.inf.ksCmH << ',' << passos1 << ','
            << e.passos - passos0 << ',' << e.newtonTotal << ',' << e.cortes << ','
            << (fim1.dt > 0. ? fim1.chuva / fim1.dt : 0.) << ',' << (fim1.dt > 0. ? fim1.dir / fim1.dt : 0.) << ','
            << (fim1.dt > 0. ? fim1.face / fim1.dt : 0.);
        std::ostringstream kv;
        kv << "caso=" << cf.caso;
        for (double th : cf.saidasH) {
            const std::string s = Num(th, 6) + "h";
            fres << ",LFS_medio_" << s << ",A_LFS_" << s << ",nivel_crista_" << s << ",frente_crista_" << s;
            const std::string a = valor(th, [](const Instante &r) { return r.m.lfsMedio; });
            const std::string b = valor(th, [](const Instante &r) { return r.m.aLFS; });
            const std::string c = valor(th, [](const Instante &r) { return r.nivelCrista; });
            const std::string d = valor(th, [](const Instante &r) { return r.frenteCrista; });
            const std::string l3 = valor(th, [](const Instante &r) { return r.m.lfsMedioLim3; });
            fres << ",LFS_medio_lim3_" << s;
            lin << ',' << a << ',' << b << ',' << c << ',' << d << ',' << l3;
            kv << " LFS_medio_" << s << '=' << a << " LFS_medio_lim3_" << s << '=' << l3 << " A_LFS_" << s << '=' << b;
        }
        const double seg = Segundos(inicio);
        fres << ",situacao,seg\n" << lin.str() << ',' << situacao << ',' << seg << '\n';
        fres.close();
        kv << " situacao=" << situacao << " seg=" << Num(seg, 6);
        {
            std::ofstream fr("resultado.txt.tmp");
            fr << kv.str() << '\n';
        }
        std::rename("resultado.txt.tmp", "resultado.txt");
        return situacao == "ok" ? 0 : 1;
    } catch (const std::exception &ex) {
        std::cerr << "V5: " << ex.what() << std::endl;
        std::ofstream fr("resultado.txt.tmp");
        fr << "situacao=erro seg=" << Num(Segundos(inicio), 6) << '\n';
        fr.close();
        std::rename("resultado.txt.tmp", "resultado.txt");
        return 1;
    }
}
