import unittest
import Correlations.choke as choke


class TestPVT(unittest.TestCase):
    def test_w_choke_gaswater_Chien_kghr(self):
        p1_atm = 20
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.w_choke_gaswater_Chien_kghr(p1_atm, p2_atm, x, gamma_g=0.55, t_C=20, d0_mm=5,
                                                                 d1_mm=100), 263, delta=10)

    def test_p_choke_up_gaswater_Chien_atm(self):
        w_kghr = 50
        p2_atm = 50
        x = 0.5
        gamma_g = 0.55
        self.assertAlmostEqual(choke.p_choke_up_gaswater_Chien_atm(w_kghr, p2_atm, x, gamma_g, t_C=20, d0_mm=5,
                                                                   d1_mm=100), 50, delta=1)

    def test_p_choke_down_gaswater_Chien_atm(self):
        w_kghr = 50
        p1_atm = 100
        x = 0.5
        gamma_g = 0.55
        self.assertAlmostEqual(choke.p_choke_down_gaswater_Chien_atm(w_kghr, p1_atm, x, gamma_g, t_C=20, d0_mm=5,
                                                                     d1_mm=100), 100, delta=1)

    def test_f_choke_difpressure_gaswater_Chien_atm(self):
        w_kghr = 500
        p1_atm = 100
        x = 0.5
        gamma_g = 0.55
        self.assertAlmostEqual(choke.f_choke_difpressure_gaswater_Chien_atm(w_kghr, p1_atm, x, gamma_g, t_C=20, d0_mm=5,
                               d1_mm=100), 3.2, delta=0.1)

    def test_w_choke_gasoil_Chien_kghr(self):
        p1_atm = 100
        p2_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.w_choke_gasoil_Chien_kghr(p1_atm, p2_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20,
                               d0_mm=5, d1_mm=100, pb_Mpa=20, rs_m3m3=300, bo_m3m3=1), 1999, delta=1)

    def test_p_choke_up_gasoil_Chien_atm(self):
        w_kghr = 50
        p2_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gasoil_Chien_atm(w_kghr, p2_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20,
                               d0_mm=5, d1_mm=100, pb_Mpa=20, rs_m3m3=300, bo_m3m3=1), 50, delta=0.1)

    def test_p_choke_down_gasoil_Chien_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gasoil_Chien_atm(w_kghr, p1_atm, x, gamma_g=0.55, gamma_oil=0.8,
                               t_C=20, d0_mm=5, d1_mm=100, pb_Mpa=20, rs_m3m3=300, bo_m3m3=1), 50,
                               delta=1)

    def test_f_choke_difpressure_gasoil_Chien_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gasoil_Chien_atm(w_kghr, p1_atm, x, gamma_g=0.55,
                               gamma_oil=0.8, t_C=20, d0_mm=5, d1_mm=100, pb_Mpa=20, rs_m3m3=300,
                               bo_m3m3=1), 0.1, delta=0.1)

    def test_w_choke_gaswater_AlSafran_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_gaswater_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C=20,
                                                                    dchoke_mm=5), 1195, delta=1)

    def test_w_choke_gasoil_AlSafran_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_gasoil_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C=20,
                                                                  dchoke_mm=5), 1216, delta=1)

    def test_w_choke_critgasoil_AlSafran_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_critgasoil_AlSafran_kghr(p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5),
                               1774, delta=1)

    def test_w_choke_critgaswater_AlSafran_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_critgaswater_AlSafran_kghr(p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5),
                               1734, delta=1)

    def test_p_choke_up_gasoil_AlSafran_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gasoil_AlSafran_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C=20,
                                                                    dchoke_mm=5), 50, delta=0.1)

    def test_p_choke_up_gaswater_AlSafran_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gaswater_AlSafran_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C=20,
                                                                      dchoke_mm=5), 50, delta=0.1)

    def test_p_choke_down_gasoil_AlSafran_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gasoil_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                                                                      dchoke_mm=5), 50, delta=0.1)

    def test_p_choke_down_gaswater_AlSafran_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gaswater_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                                                                        dchoke_mm=5), 50, delta=0.1)

    def test_f_choke_difpressure_gasoil_AlSafran_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gasoil_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                               dchoke_mm=5), 0, delta=0.1)

    def test_f_choke_difpressure_gaswater_AlSafran_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gaswater_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                                                                               dchoke_mm=5), 0, delta=0.1)

    def test_w_choke_gaswater_Perkins_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_gaswater_Perkins_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               1695, delta=1)

    def test_w_choke_gasoil_Perkins_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_gasoil_Perkins_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               1724, delta=1)

    def test_w_choke_critgaswater_Perkins_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_critgaswater_Perkins_kghr(p1_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               1695, delta=1)

    def test_w_choke_crit_gasoil_Perkins_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.w_choke_critgasoil_Perkins_kghr(p1_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               1725, delta=1)

    def test_p_choke_up_gaswater_Perkins_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gaswater_Perkins_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C=20,
                                                                     d2_mm=5), 50, delta=0.1)

    def test_p_choke_up_gasoil_Perkins_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gasoil_Perkins_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               50, delta=0.1)

    def test_p_choke_down_gaswater_Perkins_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gaswater_Perkins_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                                                                       d2_mm=5), 50, delta=0.1)

    def test_p_choke_down_gasoil_Perkins_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gasoil_Perkins_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               50, delta=0.1)

    def test_f_choke_difpressure_gaswater_Perkins_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gaswater_Perkins_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                                                                              d2_mm=5), 0, delta=0.1)

    def test_f_choke_difpressure_gasoil_Perkins_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gasoil_Perkins_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                                                                            d2_mm=5), 0, delta=0.1)

    def test_w_choke_gasoilwater_Perkins_kghr(self):
        p1_atm = 100
        p2_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.w_choke_gasoilwater_Perkins_kghr(p1_atm, p2_atm, fg, fo, fw, gamma_g=0.55, t_C=20,
                                                                      d2_mm=5, gamma_oil=0.8, pb_Mpa=20,
                                                                      rs_m3m3=300, bo_m3m3=1), 1710, delta=1)

    def test_w_choke_critgasoilwater_Perkins_kghr(self):
        p1_atm = 100
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.w_choke_critgasoilwater_Perkins_kghr(p1_atm, fg, fo, fw, gamma_g=0.55, t_C=20,
                               d2_mm=5, gamma_oil=0.8, pb_Mpa=20, rs_m3m3=300, bo_m3m3=1), 1710, delta=1)

    def test_p_choke_up_gasoilwater_Perkins_atm(self):
        w_kghr = 50
        p2_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.p_choke_up_gasoilwater_Perkins_atm(w_kghr, p2_atm, fg, fo, fw, gamma_g=0.55,
                                                                        t_C=20, d2_mm=5), 50, delta=0.1)

    def test_p_choke_down_gasoilwater_Perkins_atm(self):
        w_kghr = 50
        p1_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.p_choke_up_gasoilwater_Perkins_atm(w_kghr, p1_atm, fg, fo, fw, gamma_g=0.55,
                                                                        t_C=20, d2_mm=5), 50, delta=0.1)

    def test_f_choke_difpressure_gasoilwater_Perkins_atm(self):
        w_kghr = 50
        p1_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.f_choke_difpressure_gasoilwater_Perkins_atm(w_kghr, p1_atm, fg, fo, fw,
                               gamma_g=0.55, t_C=20, d2_mm=5), 0, delta=0.1)


    def test_w_choke_steam_Chien_kghr(self):
        p1_atm = 20
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.w_choke_steam_Chien_kghr(p1_atm, p2_atm, x,  d0_mm=5,
                                                              d1_mm=100), 225, delta=10)


    def test_p_choke_up_steam_Chien_atm(self):
        w_kghr = 50
        p2_atm = 50
        x = 1
        self.assertAlmostEqual(choke.p_choke_up_steam_Chien_atm(w_kghr, p2_atm, x, d0_mm=5, d1_mm=100),
                               50.2, delta=0.1)

    def test_p_choke_down_steam_Chien_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 1
        self.assertAlmostEqual(choke.p_choke_down_steam_Chien_atm(w_kghr, p1_atm, x, d0_mm=5, d1_mm=100),
                               49.8, delta=0.1)

    def test_f_choke_difpressure_steam_Chien_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 1
        self.assertAlmostEqual(choke.f_choke_difpressure_steam_Chien_atm(w_kghr, p1_atm, x, d0_mm=5, d1_mm=100),
                               0.2, delta=0.1)

    def test_w_choke_steam_AlSafran_kghr(self):
        p1_atm = 20
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.w_choke_steam_AlSafran_kghr(p1_atm, p2_atm, x,  dchoke_mm=5), 197, delta=10)

    def test_w_choke_critsteam_AlSafran_kghr(self):
        p1_atm = 20
        x = 1
        self.assertAlmostEqual(choke.w_choke_critsteam_AlSafran_kghr(p1_atm, x, dchoke_mm=5), 197, delta=10)

    def test_p_choke_up_steam_AlSafran_atm(self):
        w_kghr = 50
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_up_steam_AlSafran_atm(w_kghr, p2_atm, x,  dchoke_mm=5), 10.62, delta=0.1)

    def test_p_choke_down_steam_AlSafran_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_down_steam_AlSafran_atm(w_kghr, p1_atm, x,  dchoke_mm=5), 9.34, delta=0.1)

    def test_f_choke_difpressure_steam_AlSafran_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.f_choke_difpressure_steam_AlSafran_atm(w_kghr, p1_atm, x,  dchoke_mm=5),
                               0.7, delta=0.1)

    def test_w_choke_steam_Perkins_kghr(self):
        p1_atm = 20
        p2_atm = 10
        fg = 1
        self.assertAlmostEqual(choke.w_choke_steam_Perkins_kghr(p1_atm, p2_atm, fg,  d2_mm=5), 282, delta=10)

    def test_w_choke_critsteam_Perkins_kghr(self):
        p1_atm = 20
        fg = 1
        self.assertAlmostEqual(choke.w_choke_critsteam_Perkins_kghr(p1_atm, fg, d2_mm=5), 282, delta=10)

    def test_p_choke_up_steam_Perkins_atm(self):
        w_kghr = 50
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_up_steam_Perkins_atm(w_kghr, p2_atm, x,  d2_mm=5), 10.35, delta=0.1)

    def test_p_choke_down_steam_Perkins_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_down_steam_Perkins_atm(w_kghr, p1_atm, x,  d2_mm=5), 9.63, delta=0.1)

    def test_f_choke_difpressure_steamPerkins_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.f_choke_difpressure_steam_Perkins_atm(w_kghr, p1_atm, x,  d2_mm=5), 0.4,
                               delta=0.1)

if __name__ == '__main__':
    unittest.main()