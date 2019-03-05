import unittest
import Correlations.choke as choke


class TestPVT(unittest.TestCase):
    def test_w_choke_gaswater_kghr(self):
        p1_atm = 20
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.w_choke_gaswater_kghr(p1_atm, p2_atm, x, gamma_g=0.55, t_C =20, d0_mm=5,
                                d1_mm=100), 263, delta=10)

    def test_p_choke_up_gaswater_atm(self):
        w_kghr = 50
        p2_atm = 50
        x = 0.5
        gamma_g = 0.55
        self.assertAlmostEqual(choke.p_choke_up_gaswater_atm(w_kghr, p2_atm, x, gamma_g, t_C=20, d0_mm=5, d1_mm=100),
                                    50,delta=1)

    def test_p_choke_down_gaswater_atm(self):
        w_kghr = 50
        p1_atm = 100
        x = 0.5
        gamma_g = 0.55
        self.assertAlmostEqual(choke.p_choke_down_gaswater_atm(w_kghr, p1_atm, x, gamma_g, t_C=20, d0_mm=5, d1_mm=100),
                                   100,delta=1)

    def test_f_choke_difpressure_gaswater_atm(self):
        w_kghr = 500
        p1_atm = 100
        x = 0.5
        gamma_g = 0.55
        self.assertAlmostEqual(choke.f_choke_difpressure_gaswater_atm(w_kghr, p1_atm, x, gamma_g, t_C=20, d0_mm=5,
                                d1_mm=100),3.5,delta=0.1)

    def test_w_choke_gasoil_kghr(self):
        p1_atm = 100
        p2_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.w_choke_gasoil_kghr(p1_atm, p2_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20,
                                d0_mm=5, d1_mm=100, pb_Mpa=20,co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1),1975,delta=1)

    def test_p_choke_up_gasoil_atm(self):
        w_kghr = 50
        p2_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gasoil_atm(w_kghr, p2_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20,
                                d0_mm=5, d1_mm=100, pb_Mpa=20,co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1),50,delta=0.1)

    def test_p_choke_down_gasoil_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gasoil_atm(w_kghr, p1_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20,
                                    d0_mm=5, d1_mm=100,pb_Mpa=20,co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1),50,delta=1)

    def test_f_choke_difpressure_gasoil_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gasoil_atm(w_kghr, p1_atm, x, gamma_g=0.55, gamma_oil=0.8,
                            t_C=20, d0_mm=5, d1_mm=100,pb_Mpa=20,co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1),0.1,delta=0.1)

    def test_q_choke_gaswater_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.q_choke_gaswater_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C =20, dchoke_mm=5),
                               1142,delta=1)

    def test_q_choke_gasoil_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.q_choke_gasoil_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C =20, dchoke_mm=5),
                               1170,delta=1)

    def test_q_choke_critgasoil_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.q_choke_critgasoil_kghr(p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5),
                               1750, delta=1)

    def test_q_choke_critgaswater_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.q_choke_critgaswater_kghr(p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5),
                               1667, delta=1)

    def test_p_choke_up_go_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_go_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C =20, dchoke_mm=5),
                                                       50,delta=0.1)

    def test_p_choke_up_gw_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gw_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C =20, dchoke_mm=5),
                                                       50,delta=0.1)

    def test_p_choke_down_go_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_go_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5),
                               50, delta=0.1)

    def test_p_choke_down_gw_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gw_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5),
                               50, delta=0.1)

    def test_f_choke_difpressure_go_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_go_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                               dchoke_mm=5),0, delta=0.1)

    def test_f_choke_difpressure_gw_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gw_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20,
                                                                dchoke_mm=5), 0, delta=0.1)

    def test_W_choke_gaswater_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.W_choke_gaswater_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C =20, d2_mm=5),
                               1486,delta=1)

    def test_W_choke_gasoil_kghr(self):
        p1_atm = 100
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.W_choke_gasoil_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C =20, d2_mm=5),
                               1522,delta=1)

    def test_W_choke_crit_gaswater_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.W_choke_crit_gaswater_kghr(p1_atm, xg, gamma_g=0.55, t_C =20, d2_mm=5),
                               1486,delta=1)

    def test_W_choke_crit_gasoil_kghr(self):
        p1_atm = 100
        xg = 0.5
        self.assertAlmostEqual(choke.W_choke_crit_gasoil_kghr(p1_atm, xg, gamma_g=0.55, t_C =20, d2_mm=5),
                               1522,delta=1)

    def test_p_choke_up_gwater_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_gwater_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C =20, d2_mm=5),
                                                       50,delta=0.1)

    def test_p_choke_up_goil_atm(self):
        w_kghr = 50
        p2_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_up_goil_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               50, delta=0.1)

    def test_p_choke_down_gwater_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_gwater_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C =20, d2_mm=5),
                                                       50,delta=0.1)

    def test_p_choke_down_goil_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.p_choke_down_goil_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               50, delta=0.1)

    def test_f_choke_difpressure_gwater_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_gwater_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C =20, d2_mm=5),
                                                       0,delta=0.1)

    def test_f_choke_difpressure_goil_atm(self):
        w_kghr = 50
        p1_atm = 50
        xg = 0.5
        self.assertAlmostEqual(choke.f_choke_difpressure_goil_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20, d2_mm=5),
                               0, delta=0.1)

    def test_W_choke_gasoilwater_kghr(self):
        p1_atm = 100
        p2_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.W_choke_gasoilwater_kghr(p1_atm, p2_atm, fg, fo, fw, gamma_g=0.55, t_C=20, d2_mm=5,
                                gamma_oil=0.8, pb_Mpa=20, co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1), 1505, delta=1)

    def test_W_choke_crit_gasoilwater_kghr(self):
        p1_atm = 100
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.W_choke_crit_gasoilwater_kghr(p1_atm, fg, fo, fw, gamma_g=0.55, t_C=20, d2_mm=5,
                                gamma_oil=0.8, pb_Mpa=20, co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1), 1505, delta=1)

    def test_p_choke_up_gasoilwater_atm(self):
        w_kghr = 50
        p2_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.p_choke_up_gasoilwater_atm(w_kghr, p2_atm, fg, fo, fw, gamma_g=0.55, t_C =20,
                                                            d2_mm=5), 50, delta=0.1)

    def test_p_choke_down_gasoilwater_atm(self):
        w_kghr = 50
        p1_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.p_choke_up_gasoilwater_atm(w_kghr, p1_atm, fg, fo, fw, gamma_g=0.55, t_C =20,
                                                                d2_mm=5), 50, delta=0.1)

    def test_f_choke_difpressure_gasoilwater_atm(self):
        w_kghr = 50
        p1_atm = 50
        fg = 0.5
        fo = 0.3
        fw = 0.2
        self.assertAlmostEqual(choke.f_choke_difpressure_gasoilwater_atm(w_kghr, p1_atm, fg, fo, fw, gamma_g=0.55,
                                                t_C =20, d2_mm=5), 0, delta=0.1)


    def test_w_choke_saturatedsteam_kghr(self):
        p1_atm = 20
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.w_choke_saturatedsteam_kghr(p1_atm, p2_atm, x,  d0_mm=5,
                                d1_mm=100), 225, delta=10)


    def test_p_choke_up_saturatedsteam_atm(self):
        w_kghr = 50
        p2_atm = 50
        x = 1
        self.assertAlmostEqual(choke.p_choke_up_saturatedsteam_atm(w_kghr, p2_atm, x, d0_mm=5, d1_mm=100),
                               50.2, delta=0.1)

    def test_p_choke_down_saturatedsteam_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 1
        self.assertAlmostEqual(choke.p_choke_down_saturatedsteam_atm(w_kghr, p1_atm, x, d0_mm=5, d1_mm=100),
                               49.8, delta=0.1)

    def test_f_choke_difpressure_saturatedsteam_atm(self):
        w_kghr = 50
        p1_atm = 50
        x = 1
        self.assertAlmostEqual(choke.f_choke_difpressure_saturatedsteam_atm(w_kghr, p1_atm, x, d0_mm=5, d1_mm=100),
                               0.2, delta=0.1)

    def test_q_choke_saturatsteam_kghr(self):
        p1_atm = 20
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.q_choke_saturatsteam_kghr(p1_atm, p2_atm, x,  dchoke_mm=5), 197, delta=10)

    def test_q_choke_critsaturatsteam_kghr(self):
        p1_atm = 20
        x = 1
        self.assertAlmostEqual(choke.q_choke_critsaturatsteam_kghr(p1_atm, x, dchoke_mm=5), 197, delta=10)

    def test_p_choke_up_saturatsteam_atm(self):
        w_kghr = 50
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_up_saturatsteam_atm(w_kghr, p2_atm, x,  dchoke_mm=5), 10.62, delta=0.1)

    def test_p_choke_down_saturatsteam_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_down_saturatsteam_atm(w_kghr, p1_atm, x,  dchoke_mm=5), 9.34, delta=0.1)

    def test_f_choke_difpressure_saturatsteam_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.f_choke_difpressure_saturatsteam_atm(w_kghr, p1_atm, x,  dchoke_mm=5),
                               0.7, delta=0.1)

    def test_W_choke_ssteam_kghr(self):
        p1_atm = 20
        p2_atm = 10
        fg = 1
        self.assertAlmostEqual(choke.W_choke_ssteam_kghr(p1_atm, p2_atm, fg,  d2_mm=5), 282, delta=10)

    def test_W_choke_critssteam_kghr(self):
        p1_atm = 20
        fg = 1
        self.assertAlmostEqual(choke.W_choke_critssteam_kghr(p1_atm, fg, d2_mm=5), 282, delta=10)

    def test_p_choke_up_ssteam_atm(self):
        w_kghr = 50
        p2_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_up_ssteam_atm(w_kghr, p2_atm, x,  d2_mm=5), 10.35, delta=0.1)

    def test_p_choke_down_ssteam_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.p_choke_down_ssteam_atm(w_kghr, p1_atm, x,  d2_mm=5), 9.63, delta=0.1)

    def test_f_choke_difpressure_ssteam_atm(self):
        w_kghr = 50
        p1_atm = 10
        x = 1
        self.assertAlmostEqual(choke.f_choke_difpressure_ssteam_atm(w_kghr, p1_atm, x,  d2_mm=5), 0.4, delta=0.1)

if __name__ == '__main__':
    unittest.main()