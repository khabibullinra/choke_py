"""
файл для отладки расчета штуцера с газом и водой
"""

import Correlations.choke as choke


print(choke.w_choke_gaswater_kghr(p1_atm=10, p2_atm=2, x=0.1, gamma_g=0.55, t_C =20, d0_mm=5, d1_mm=100))

#тут можно сделать, чтобы рисовались графики показывающие, что считается
