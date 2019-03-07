"""
файл для отладки расчета штуцера с паром
"""
import numpy as np
import matplotlib.pyplot as plt
import Correlations.choke as choke

def plot_w_choke_steam(p1_atm):
    p2_a = np.arange(0.1, p1_atm, 0.05)
    q1 = np.array([])
    q2 = np.array([])
    q3 = np.array([])
    for p2_atm in p2_a:
        p1 = choke.w_choke_steam_Chien_kghr(p1_atm, p2_atm, 1, d0_mm=5, d1_mm=100)
        p2 = choke.w_choke_steam_AlSafran_kghr(p1_atm, p2_atm, 1, dchoke_mm=5)
        p3 = choke.w_choke_steam_Perkins_kghr(p1_atm, p2_atm, 1, d2_mm=5)
        q1 = np.append(q1, p1)
        q2 = np.append(q2, p2)
        q3 = np.append(q3, p3)
    plt.plot(p2_a, q3, label='Перкинс')
    plt.plot(p2_a, q2, label='Альсафран')
    plt.plot(p2_a, q1, label='Чиен (Миллер)')
    plt.title('Расчет расхода пара от давления на выходе разными методиками')
    plt.ylabel('Расход пара, кг/час')
    plt.xlabel('Давление на выходе, атм')
    plt.legend()
    return plt


'Построение графика "Давление на входе в зависимости от расхода и давления на выходе"'
def plot_p_choke_up_steam(p2_atm, w_kghr_0, w_kghr_max):
    w_kghr_a = np.arange(w_kghr_0, w_kghr_max, 2)  # расход, кг/час для которого находится давление на выходе
    q1 = np.array([])
    q2 = np.array([])
    q3 = np.array([])
    for w_kghr in w_kghr_a:
        p1 = choke.p_choke_up_steam_Chien_atm(w_kghr, p2_atm, 1)
        p2 = choke.p_choke_up_steam_AlSafran_atm(w_kghr, p2_atm, 1)
        p3 = choke.p_choke_up_steam_Perkins_atm(w_kghr, p2_atm, 1)
        q1 = np.append(q1, p1)
        q2 = np.append(q2, p2)
        q3 = np.append(q3, p3)
    plt.plot(w_kghr_a, q3, label='Перкинс')
    plt.plot(w_kghr_a, q2, label='Альсафран')
    plt.plot(w_kghr_a, q1, label='Чиен (Миллер)')
    plt.title('Давление на входе в зависимости от расхода и давления на выходе')
    plt.xlabel('Расход пара, кг/час')
    plt.ylabel('Давление на входе, атм')
    plt.legend()
    return plt


'Построение графика "Давление на выходе в зависимости от расхода и давления на входе"'
def plot_choke_down_steam(p1_atm, w_kghr_0, w_kghr_max):
    w_kghr_a = np.arange(w_kghr_0, w_kghr_max, 2)  # расход, кг/час для которого находится давление на выходе
    q1 = np.array([])
    q2 = np.array([])
    q3 = np.array([])
    for w_kghr in w_kghr_a:
        p1 = choke.p_choke_down_steam_Chien_atm(w_kghr, p1_atm, 1)
        p2 = choke.p_choke_down_steam_AlSafran_atm(w_kghr, p1_atm, 1)
        p3 = choke.p_choke_down_steam_Perkins_atm(w_kghr, p1_atm, 1)
        q1 = np.append(q1, p1)
        q2 = np.append(q2, p2)
        q3 = np.append(q3, p3)
    plt.plot(w_kghr_a, q3, label='Перкинс')
    plt.plot(w_kghr_a, q2, label='Альсафран')
    plt.plot(w_kghr_a, q1, label='Чиен (Миллер)')
    plt.title('Давление на выходе в зависимости от расхода и давления на входе')
    plt.xlabel('Расход пара, кг/час')
    plt.ylabel('Давление на выходе, атм')
    plt.legend()
    return plt
