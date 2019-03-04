"""
файл для отладки расчета штуцера с нефтью, газом и водой
"""
import numpy as np
import matplotlib.pyplot as plt
import Correlations.choke as choke



p1_atm = 20
p2_a = np.arange(0.1, p1_atm, 0.05)
q3 = np.array([])
for p2_atm in p2_a:
    p3 = choke.W_choke_gasoilwater_kghr(p1_atm, p2_atm, 1, 0, 0)
    q3 = np.append(q3, p3)
plt.plot(p2_a, q3, label='Перкинс')
plt.title('Расчет расхода водонефтегазовой смеси от давления на выходе')
plt.ylabel('Расход смеси, кг/час')
plt.xlabel('Давление на выходе, атм')
plt.legend()
plt.show()

'Построение графика "Давление на входе в зависимости от расхода и давления на выходе"'
p2_atm = 10
w_kghr_a = np.arange(0, 200, 2) # расход, кг/час для которого находится давление на выходе
q3 = np.array([])
for w_kghr in w_kghr_a:
    p3 = choke.p_choke_up_gasoilwater_atm(w_kghr, p2_atm, 1,0,0)
    q3 = np.append(q3, p3)
plt.plot(w_kghr_a, q3, label='Перкинс')
plt.title('Давление на входе в зависимости от расхода и давления на выходе')
plt.xlabel('Расход смеси, кг/час')
plt.ylabel('Давление на входе, атм')
plt.legend()
plt.show()


'Построение графика "Давление на выходе в зависимости от расхода и давления на входе"'
p1_atm = 20
w_kghr_a = np.arange(0, 200, 2) # расход, кг/час для которого находится давление на выходе
q3 = np.array([])
for w_kghr in w_kghr_a:
    p3 = choke.p_choke_down_gasoilwater_atm(w_kghr, p1_atm, 1,0,0)
    q3 = np.append(q3, p3)
plt.plot(w_kghr_a, q3, label='Перкинс')
plt.title('Давление на выходе в зависимости от расхода и давления на входе')
plt.xlabel('Расход смеси, кг/час')
plt.ylabel('Давление на выходе, атм')
plt.legend()
plt.show()