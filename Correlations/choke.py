import scipy.optimize as opt
import uPVT.PVT as PVT
from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)


def w_choke_gaswater_Chien_kghr(p1_atm, p2_atm, x, gamma_g=0.55, t_C=20, d0_mm=5, d1_mm=100):
    """
    расчет расхода водогазовой смеси через штуцер по методике Чиена (Миллера)
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3

    c0 = 0.8  # коэффициент разряда штуцера
    y1 = 0.89487  # коэффициент расширения пара,считают его постоянным,зависит от давлений и диаметров
    a = 0.99998  # эмпирические коэффициенты, рассчитаны из графиков по статье
    b = 1.38
    vfg_m3kg = 1 / rog_kgm3  # удельный объем газа при заданном давлении м3/кг
    vf_m3kg = steamTable.vL_p(p1_atm)  # удельный объем воды при заданном давлении м3/кг
    vexp_m3kg = a * vfg_m3kg * x ** b + vf_m3kg  # удельный объем смеси м3/кг
    ro_kgm3 = 1 / vexp_m3kg  # плотность смеси кг/м3
    beta = d0_mm / d1_mm
    f_atm = p1_atm - p2_atm
    f_kpa = f_atm * 100
    w_kgsec = 3.512407 * 10 ** (-5) * c0 * y1 * d0_mm ** 2 * (f_kpa * ro_kgm3) ** 0.5 / (1 - beta ** 4) ** 0.5
    w_kghr = w_kgsec * 3600
    return w_kghr



def p_choke_up_gaswater_Chien_atm(w_kghr, p2_atm, x, gamma_g, t_C=20, d0_mm=5, d1_mm=100):
    """

    расчет давления водогазовой смеси перед клапаном по методике Чиена (Миллера)
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    """

    def w2(p1_atm):
        return w_choke_gaswater_Chien_kghr(p1_atm, p2_atm, x, gamma_g, t_C, d0_mm, d1_mm) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_gaswater_Chien_atm(w_kghr, p1_atm, x, gamma_g, t_C=20, d0_mm=5, d1_mm=100):
    """
    расчет давления на выходе от давления на входе при разных расходах по методике Чиена (Миллера)
    p1_atm - давление на входе в штуцер, атм
    w_kghr- расход пара, кг/час
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    d0_mm - диаметр штуцера, мм
    d1_mm - диаметр трубы, мм
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3

    c0 = 0.8  # коэффициент разряда штуцера
    y1 = 0.89487  # коэффициент расширения пара
    a = 0.99998  # эмпирические коэффициенты, рассчитаны из графиков по статье
    b = 1.38
    vfg_m3kg = 1 / rog_kgm3  # удельный объем газа при заданном давлении м3/кг
    vf_m3kg = steamTable.vL_p(p1_atm)  # удельный объем воды при заданном давлении м3/кг
    vexp_m3kg = a * vfg_m3kg * x ** b + vf_m3kg  # удельный объем смеси м3/кг
    ro_kgm3 = 1 / vexp_m3kg  # плотность смеси кг/м3
    beta = d0_mm / d1_mm
    w_kgsec = w_kghr / 3600
    f_kpa = w_kgsec ** 2 * (1 - beta ** 4) * 10 ** 10 / (12.337 * c0 ** 2 * y1 ** 2 * d0_mm ** 4 * ro_kgm3)
    p2_atm = p1_atm - f_kpa / 100
    return p2_atm


def f_choke_difpressure_gaswater_Chien_atm(w_kghr, p1_atm, x, gamma_g, t_C=20, d0_mm=5, d1_mm=100):
    """
    расчет разницы давлений до и после штуцера водогазовой смеси по методике Чиена (Миллера)
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    """

    def w2(f):
        return w_choke_gaswater_Chien_kghr(p1_atm, p1_atm - f, x, gamma_g, t_C, d0_mm, d1_mm) - w_kghr

    return opt.fsolve(w2, 0)


def w_choke_gasoil_Chien_kghr(p1_atm, p2_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20, d0_mm=5, d1_mm=100, pb_Mpa=20,
                              rs_m3m3=300, bo_m3m3=1):
    """
    расчет расхода нефтегазовой смеси через штуцер по методике Чиена (Миллера)
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    gamma_oil - относительная плонтость нефти по воде
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3

    c0 = 0.8  # коэффициент разряда штуцера
    y1 = 0.89487  # коэффициент расширения пара,считают его постоянным,зависит от давлений и диаметров
    a = 0.99998  # эмпирические коэффициенты, рассчитаны из графиков по статье
    b = 1.38
    vfg_m3kg = 1 / rog_kgm3  # удельный объем газа при заданном давлении м3/кг
    vf_m3kg = 1 / roo_kgm3  # удельный объем нефти при заданном давлении м3/кг
    vexp_m3kg = a * vfg_m3kg * x ** b + vf_m3kg  # удельный объем смеси м3/кг
    ro_kgm3 = 1 / vexp_m3kg  # плотность смеси кг/м3
    beta = d0_mm / d1_mm
    f_atm = p1_atm - p2_atm
    f_kpa = f_atm * 100
    w_kgsec = 3.512407 * 10 ** (-5) * c0 * y1 * d0_mm ** 2 * (f_kpa * ro_kgm3) ** 0.5 / (1 - beta ** 4) ** 0.5
    w_kghr = w_kgsec * 3600
    return w_kghr



def p_choke_up_gasoil_Chien_atm (w_kghr, p2_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20, d0_mm=5, d1_mm=100, pb_Mpa=20,
                                 rs_m3m3=300, bo_m3m3=1):
    """

    расчет давления водогазовой смеси перед клапаном по методике Чиена (Миллера)
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    gamma_oil - относительная плонтость нефти по воде
    x - массовая доля газа, д.ед.
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(p1_atm):
        return w_choke_gasoil_Chien_kghr(p1_atm, p2_atm, x, gamma_g, gamma_oil, t_C, d0_mm, d1_mm, pb_Mpa,
                                         rs_m3m3, bo_m3m3) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_gasoil_Chien_atm(w_kghr, p1_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20, d0_mm=5, d1_mm=100,
                                  pb_Mpa=20, rs_m3m3=300, bo_m3m3=1):
    """
    расчет давления на выходе от давления на входе при разных расходах по методике Чиена (Миллера)
    p1_atm - давление на входе в штуцер, атм
    w_kghr- расход пара, кг/час
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    d0_mm - диаметр штуцера, мм
    d1_mm - диаметр трубы, мм
    gamma_oil - относительная плонтость нефти по воде
    x - массовая доля газа, д.ед.
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """
    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3

    c0 = 0.8  # коэффициент разряда штуцера
    y1 = 0.89487  # коэффициент расширения пара
    a = 0.99998  # эмпирические коэффициенты, рассчитаны из графиков по статье
    b = 1.38
    vfg_m3kg = 1 / rog_kgm3  # удельный объем газа при заданном давлении м3/кг
    vf_m3kg = 1 / roo_kgm3  # удельный объем нефти при заданном давлении м3/кг
    vexp_m3kg = a * vfg_m3kg * x ** b + vf_m3kg  # удельный объем смеси м3/кг
    ro_kgm3 = 1 / vexp_m3kg  # плотность смеси кг/м3
    beta = d0_mm / d1_mm
    w_kgsec = w_kghr / 3600
    f_kpa = w_kgsec ** 2 * (1 - beta ** 4) * 10 ** 10 / (12.337 * c0 ** 2 * y1 ** 2 * d0_mm ** 4 * ro_kgm3)
    p2_atm = p1_atm - f_kpa / 100
    return p2_atm


def f_choke_difpressure_gasoil_Chien_atm(w_kghr, p1_atm, x, gamma_g=0.55, gamma_oil=0.8, t_C=20, d0_mm=5, d1_mm=100, 
                                         pb_Mpa=20, rs_m3m3=300, bo_m3m3=1):
    """
    расчет разницы давлений до и после штуцера водогазовой смеси по методике Чиена (Миллера)
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    gamma_oil - относительная плонтость нефти по воде
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(f):
        return w_choke_gasoil_Chien_kghr(p1_atm, p1_atm - f, x, gamma_g, gamma_oil, t_C, d0_mm, d1_mm, pb_Mpa,
                                         rs_m3m3, bo_m3m3) - w_kghr

    return opt.fsolve(w2, 0)


def w_choke_gaswater_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5):
    """
    расчет расхода водогазовой смеси через штуцер по методике Альсафран
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    dchoke_mm - диаметр штуцера,мм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    """
    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3


    if xg <= 0:xg = 0.0001
    cd = 0.75  # дисчарч коэффициент для штуцера (коэффициент разряда штуцера)
    g_c = 9.8  # ускорение свободного падения
    p1_pa = p1_atm * 100  # давление до штуцера, Па
    p2_pa = p2_atm * 100  # давление после штуцера, Па
    vgl = 1 / rog_kgm3
    vl = steamTable.vL_p(p1_atm)  # удельный объем воды, м3/кг
    cl = steamTable.CvL_p(p1_atm)  # удельная теплоемкость воды, кДж/кг/К
    cvg = PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    k = cpg / cvg  # показатель адиабаты
    d2 = dchoke_mm / 1000  # диаметр клапана, м
    r = (1 + xg * (vgl / vl - 1)) ** 0.5 * (1 + 0.6 / 2.718 ** (5 * xg))
    alfa = r * (1 - xg) * vl / (xg * vgl)
    a2 = 3.14 * d2 ** 2 / 4
    n = (xg * k * cvg + (1 - xg) * cl) / (xg * cvg + (1 - xg) * cl)

    def y(x):
        y = (alfa * (1 - x) + n / (n - 1)) / (n / (n - 1) + n / 2 * (1 + alfa * x ** (1 / n)) ** 2) - x ** (1 - 1 / n)
        return (y)

    rc = opt.fsolve(y, 0.5)
    r= p2_pa / p1_pa
    if rc >= r:
        pr = rc
    else:
        pr = r
    if pr > 1: pr = 1
    ab = n / (n - 1) * (1 - pr ** ((n - 1) / n)) + alfa * (1 - pr)
    ac = xg * vgl * (pr ** (-1 / n) + alfa) ** 2 * (xg + 1 / r * (1 - xg))
    wi = a2 * (288 * g_c * cd ** 2 * p1_pa * ab / ac) ** 0.5 * 3600
    return wi


def w_choke_gasoil_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g=0.55, gamma_oil=0.8, t_C=20, dchoke_mm=5, pb_Mpa=20,
                                 co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет расхода нефтегазовой смеси через штуцер по методике Альсафран
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    dchoke_mm - диаметр штуцера,мм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3


    if xg <= 0: xg = 0.0001
    cd = 0.75  # дисчарч коэффициент для штуцера (коэффициент разряда штуцера)
    g_c = 9.8  # ускорение свободного падения
    p1_pa = p1_atm * 100  # давление до штуцера, Па
    p2_pa = p2_atm * 100  # давление после штуцера, Па
    vgl = 1 / rog_kgm3
    vl = 1 / roo_kgm3
    cl = PVT.Cvo_kJkgK  # удельная теплоемкость нефти, кДж/кг/К
    cvg = PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    k = cpg / cvg # показатель адиабаты
    d2 = dchoke_mm / 1000  # диаметр клапана, м
    r = (1 + xg * (vgl / vl - 1)) ** 0.5 * (1 + 0.6 / 2.718 ** (5 * xg))
    alfa = r * (1 - xg) * vl / (xg * vgl)
    a2 = 3.14 * d2 ** 2 / 4
    n = (xg * k * cvg + (1 - xg) * cl) / (xg * cvg + (1 - xg) * cl)

    def y(x):
        y = (alfa * (1 - x) + n / (n - 1)) / (n / (n - 1) + n / 2 * (1 + alfa * x ** (1 / n)) ** 2) - x ** (1 - 1 / n)
        return (y)

    rc = opt.fsolve(y, 0.5)
    r = p2_pa / p1_pa
    if rc >= r:
        pr = rc
    else:
        pr = r
    if pr > 1: pr = 1
    ab = n / (n - 1) * (1 - pr ** ((n - 1) / n)) + alfa * (1 - pr)
    ac = xg * vgl * (pr ** (-1 / n) + alfa) ** 2 * (xg + 1 / r * (1 - xg))
    wi = a2 * (288 * g_c * cd ** 2 * p1_pa * ab / ac) ** 0.5 * 3600
    return wi


def w_choke_critgaswater_AlSafran_kghr(p1_atm, xg, dchoke_mm=5, gamma_g=0.55, t_C=20):
    """
    расчет критического расхода водогазовой смеси через штуцер по методике Альсафран
    p1_atm - давление на входе в штуцер, атм
    dchoke_mm - диаметр штуцера,мм
    xg - массовая доля газа, д.ед.
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3


    if xg <= 0: xg = 0.0001
    cd = 0.75  # дисчарч коэффициент для штуцера (коэффициент разряда штуцера)
    g_c = 9.8  # ускорение свободного падения
    p1_pa = p1_atm * 100  # давление до штуцера, Па
    vgl = 1 / rog_kgm3
    vl = steamTable.vL_p(p1_atm)  # удельный объем воды, м3/кг
    cl = steamTable.CvL_p(p1_atm)  # удельная теплоемкость воды, кДж/кг/К
    cvg = PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    k = cpg / cvg # показатель адиабаты
    d2 = dchoke_mm / 1000  # диаметр клапана, м
    r = (1 + xg * (vgl / vl - 1)) ** 0.5 * (1 + 0.6 / 2.718 ** (5 * xg))
    alfa = r * (1 - xg) * vl / (xg * vgl)
    a2 = 3.14 * d2 ** 2 / 4
    n = (xg * k * cvg + (1 - xg) * cl) / (xg * cvg + (1 - xg) * cl)

    def y(x):
        y = (alfa * (1 - x) + n / (n - 1)) / (n / (n - 1) + n / 2 * (1 + alfa * x ** (1 / n)) ** 2) - x ** (1 - 1 / n)
        return (y)

    rc = opt.fsolve(y, 0.5)
    pr = rc
    if pr > 1:
        pr = 1
    ab = n / (n - 1) * (1 - pr ** ((n - 1) / n)) + alfa * (1 - pr)
    ac = xg * vgl * (pr ** (-1 / n) + alfa) ** 2 * (xg + 1 / r * (1 - xg))
    wi = a2 * (288 * g_c * cd ** 2 * p1_pa * ab / ac) ** 0.5 * 3600
    return wi


def w_choke_critgasoil_AlSafran_kghr(p1_atm, xg, gamma_g=0.55, gamma_oil=0.8, t_C =20, dchoke_mm=5, pb_Mpa=20, 
                                     co_1Mpa =0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет критического расхода нефтегазовой смеси через штуцер по методике Альсафран
    p1_atm - давление на входе в штуцер, атм
    dchoke_mm - диаметр штуцера,мм
    xg - массовая доля газа, д.ед.
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """
    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3

    if xg <= 0: xg = 0.0001
    cd = 0.75  # дисчарч коэффициент для штуцера (коэффициент разряда штуцера)
    g_c = 9.8  # ускорение свободного падения
    p1_pa = p1_atm * 100  # давление до штуцера, Па
    vgl = 1 / rog_kgm3  # удельный объем газа, м3/кг
    vl = 1 / roo_kgm3
    cl = PVT.Cvo_kJkgK  # удельная теплоемкость нефти, кДж/кг/К
    cvg = PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    k = cpg / cvg   # показатель адиабаты
    d2 = dchoke_mm / 1000  # диаметр клапана, м
    r = (1 + xg * (vgl / vl - 1)) ** 0.5 * (1 + 0.6 / 2.718 ** (5 * xg))
    alfa = r * (1 - xg) * vl / (xg * vgl)
    a2 = 3.14 * d2 ** 2 / 4
    n = (xg * k * cvg + (1 - xg) * cl) / (xg * cvg + (1 - xg) * cl)

    def y(x):
        y = (alfa * (1 - x) + n / (n - 1)) / (n / (n - 1) + n / 2 * (1 + alfa * x ** (1 / n)) ** 2) - x ** (1 - 1 / n)
        return (y)

    rc = opt.fsolve(y, 0.5)
    pr = rc
    if pr > 1: pr = 1
    ab = n / (n - 1) * (1 - pr ** ((n - 1) / n)) + alfa * (1 - pr)
    ac = xg * vgl * (pr ** (-1 / n) + alfa) ** 2 * (xg + 1 / r * (1 - xg))
    wi = a2 * (288 * g_c * cd ** 2 * p1_pa * ab / ac) ** 0.5 * 3600
    return wi


def p_choke_up_gaswater_AlSafran_atm(w_kghr, p2_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5):
    """
    расчет давления водогазовой смеси перед клапаном по методике Альсафран
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    dchoke_mm - диаметр штуцера, мм
    """

    def w2(p1_atm):
        return w_choke_gaswater_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g, t_C, dchoke_mm) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_up_gasoil_AlSafran_atm(w_kghr, p2_atm, xg, gamma_g=0.55, gamma_oil=0.8, t_C=20, dchoke_mm=5, pb_Mpa=20,
                                   co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет давления нефтегазовой смеси перед клапаном по методике Альсафран
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    dchoke_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(p1_atm):
        return w_choke_gasoil_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g, gamma_oil, t_C, dchoke_mm, pb_Mpa,
                                            co_1Mpa, rs_m3m3, bo_m3m3) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_gaswater_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5):
    """
    расчет давления водогазовой смеси перед клапаном по методике Альсафран
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    dchoke_mm - диаметр штуцера, мм
    """

    def w2(p2_atm):
        return w_choke_gaswater_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g, t_C, dchoke_mm) - w_kghr

    p2 = opt.fsolve(w2, p1_atm - 0.00001)
    if w_kghr > w_choke_critgaswater_AlSafran_kghr(p1_atm, xg, dchoke_mm, gamma_g, t_C):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        p2 = 0
    return p2


def p_choke_down_gasoil_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, gamma_oil=0.8, t_C=20, dchoke_mm=5, pb_Mpa=20,
                                     co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет давления газонефтяной смеси перед клапаном по методике Альсафран
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    dchoke_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(p2_atm):
        return w_choke_gasoil_AlSafran_kghr(p1_atm, p2_atm, xg, gamma_g, gamma_oil, t_C, dchoke_mm, pb_Mpa,
                                            co_1Mpa, rs_m3m3, bo_m3m3) - w_kghr

    p2 = opt.fsolve(w2, p1_atm - 0.00001)
    if w_kghr > w_choke_critgasoil_AlSafran_kghr(p1_atm, xg, gamma_g, gamma_oil, t_C, dchoke_mm, pb_Mpa,
                                                 co_1Mpa, rs_m3m3, bo_m3m3):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        p2 = 0
    return p2


def f_choke_difpressure_gaswater_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, t_C=20, dchoke_mm=5):
    """
    расчет разницы давлений до и после штуцера водогазовой смеси по методике Aльсафран
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    dchoke_mm - диаметр штуцера, мм
    """

    def w2(f):
        return w_choke_gaswater_AlSafran_kghr(p1_atm, p1_atm - f, xg, gamma_g, t_C, dchoke_mm) - w_kghr

    f = opt.fsolve(w2, 0)
    if w_kghr > w_choke_critgaswater_AlSafran_kghr(p1_atm, xg, dchoke_mm, gamma_g, t_C):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        f = 0
    return f


def f_choke_difpressure_gasoil_AlSafran_atm(w_kghr, p1_atm, xg, gamma_g=0.55, gamma_oil=0.8, t_C=20, dchoke_mm=5, 
                                            pb_Mpa=20, co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет разницы давлений до и после штуцера нефтегазовой смеси по методике Aльсафран
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    dchoke_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(f):
        return w_choke_gasoil_AlSafran_kghr(p1_atm, p1_atm - f, xg, gamma_g, gamma_oil, t_C, dchoke_mm, pb_Mpa,
                                            co_1Mpa, rs_m3m3, bo_m3m3) - w_kghr

    f = opt.fsolve(w2, 0)
    if w_kghr > w_choke_critgasoil_AlSafran_kghr(p1_atm, xg, gamma_g, gamma_oil, t_C, dchoke_mm, pb_Mpa,
                                                 co_1Mpa, rs_m3m3, bo_m3m3):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        f = 0
    return f


def w_choke_gaswater_Perkins_kghr (p1_atm, p2_atm, fg, gamma_g=0.55, t_C=20, d2_mm=5):
    """
    расчет расхода водогазовой смеси через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    t_C - температура, С
    d2_mm - диаметр штуцера, мм
    gamma_g - относительная плотность газа по воздуху
    fg - весовая доля газа в потоке, д.ед.
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3

    if fg <= 0: fg = 0.00001
    fw = 1 - fg  # весовая доля воды в потоке
    m_m = 16  # молекулярный вес пара, моль газа
    z = fl.z  # коэффициент сжимаемости газа
    row_kgm3 = steamTable.rhoL_p(p1_atm)  # плотность воды   кг/м3
    cvw = 0.24 * 778.169 * steamTable.CvL_p(p1_atm)  # удельная теплоемкость воды 
    # 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
    cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    row = row_kgm3 * 0.062428  # плотность воды  Ibm/ft3
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    f = cpg / cvg   # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    p2_psia = p2_atm * 14.2233  # давление за штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * fw / row
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fw * cvw) / (fg * cvg + fw * cvw)
    lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
    i = 0
    pr1 = 0.1
    pr11 = 0.2
    """
    Пересчитываю отношение давлений столько раз, пока разница между следующим и предыдущим не будет очень маленькой.
    Отношение давлений служит для того, чтобы найти давление в штуцере.
    """
    while abs(pr1 - pr11) > 0.001 and i < 10:  # пересчитываем удельные теплоемкости и давление в штуцере, пока
                                                                    #  они не будут соответствовать A*(B+C)=D*E)
        
        def qw(pr2):
            a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
            b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        fg / n * pr2 ** (-(1 + n) / n))
            c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                        (fg * pr2 ** (-1 / n)) + alf) ** 2
            d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        (fg * pr2 ** (-1 / n)) + alf)
            e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
            return a * (b + c) - d * e

        pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
        p3_psia = p1_psia * pr1  # находим давление в штуцере
        pmid_atm = (p1_psia + p3_psia) / 2 / 14.2233  # среднее давление в штуцере и перед штуцером, пересчитываем 
        # удельные теплоемкости для этого давления
        cvw = 0.24 * 778.169 * steamTable.CvL_p(pmid_atm)  # удельная теплоемкость воды 
        # 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
        cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
        cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
        f = cpg / cvg
        n = (fg * f * cvg + fw * cvw) / (fg * cvg + fw * cvw)
        lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr ))
        pr11 = opt.fsolve(qw, 0.01)
        i = i + 1
    if p3_psia > p2_psia:
        pr = p3_psia / p1_psia
    else:
        pr = p2_psia / p1_psia
    if pr > 1: pr = 1
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
            (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) *
                        ((fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia / v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def w_choke_gasoil_Perkins_kghr(p1_atm, p2_atm, fg, gamma_g=0.55, t_C=20, d2_mm=5, gamma_oil=0.8, pb_Mpa =20,
                                co_1Mpa =0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет расхода нефтегазовой смеси через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    t_C - температура, С
    d2_mm - диаметр штуцера, мм
    gamma_g - относительная плотность газа по воздуху
    fg - весовая доля газа в потоке, д.ед.
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3

    if fg <= 0: fg = 0.00001
    fo = 1 - fg  # весовая доля нефти в потоке
    m_m = 16  # молекулярный вес пара, моль газа
    z = fl.z  # коэффициент сжимаемости газа
    cvo = 0.24 * 778.169 * PVT.Cvo_kJkgK  # удельная теплоемкость нефти, кДж/кг/К
    cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    roo = roo_kgm3 * 0.062428  # плотность воды  Ibm/ft3
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    f = cpg / cvg   # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    p2_psia = p2_atm * 14.2233  # давление за штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * fo / roo
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fo * cvo) / (fg * cvg + fo * cvo)
    lambd = fg + ((fg * cvg + fo * cvo) * m_m / (z * r_ftIbtlbmmolr))

    def qw(pr2):
        a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
        b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                    fg / n * pr2 ** (-(1 + n) / n))
        c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                    (fg * pr2 ** (-1 / n)) + alf) ** 2
        d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                (fg * pr2 ** (-1 / n)) + alf)
        e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
        return a * (b + c) - d * e

    pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
    p3_psia = p1_psia * pr1  # находим давление в штуцере
    if p3_psia > p2_psia:
        pr = p3_psia / p1_psia
    else:
        pr = p2_psia / p1_psia
    if pr > 1: pr = 1
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
                (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia / v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def w_choke_critgaswater_Perkins_kghr(p1_atm, fg, gamma_g=0.55, t_C=20, d2_mm=5):
    """
    расчет критического расхода водогазовой смеси через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    t_C - температура, С
    gamma_g - относительная плотность газа по воздуху
    d2_mm - диаметр штуцера, мм
    fg - весовая доля газа в потоке, д.ед.
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3


    if fg <= 0: fg = 0.00001
    fw = 1 - fg  # весовая доля воды в потоке
    m_m = 16  # молекулярный вес пара, моль газа
    z = fl.z  # коэффициент сжимаемости газа
    row_kgm3 = steamTable.rhoL_p(p1_atm)  # плотность воды   кг/м3
    cvw = 0.24 * 778.169 * steamTable.CvL_p(
        p1_atm)  # удельная теплоемкость воды 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
    cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    row = row_kgm3 * 0.062428  # плотность воды  Ibm/ft3
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    f = cpg / cvg # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * fw / row
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fw * cvw) / (fg * cvg + fw * cvw)
    lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
    i = 0
    pr1 = 0.1
    pr11 = 0.2
    """
    Пересчитываю отношение давлений столько раз, пока разница между следующим и предыдущим не будет очень маленькой.
    Отношение давлений служит для того, чтобы найти давление в штуцере.
    """
    while abs(pr1 - pr11) > 0.001 and i < 10:  # пересчитываем удельные теплоемкости и давление в штуцере, 
        # пока они не будут соответствовать A*(B+C)=D*E
        def qw(pr2):
            a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
            b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        fg / n * pr2 ** (-(1 + n) / n))
            c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                    (fg * pr2 ** (-1 / n)) + alf) ** 2
            d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        (fg * pr2 ** (-1 / n)) + alf)
            e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
            return a * (b + c) - d * e

        pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
        p3_psia = p1_psia * pr1  # находим давление в штуцере
        pmid_atm = (p1_psia + p3_psia) / 2 / 14.2233  # среднее давление в штуцере и перед штуцером, 
    # пересчитываем удельные теплоемкости для этого давления
        cvw = 0.24 * 778.169 * steamTable.CvL_p(
            pmid_atm)  # удельная теплоемкость жидкости 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
        cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
        cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
        f = cpg / cvg
        n = (fg * f * cvg + fw * cvw) / (fg * cvg + fw * cvw)
        lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr ))
        pr11 = opt.fsolve(qw, 0.01)
        i = i + 1
    pr = p3_psia / p1_psia
    if pr > 1: pr = 1
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
                (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia / v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def w_choke_critgasoil_Perkins_kghr(p1_atm, fg, gamma_g=0.55, t_C=20, d2_mm=5, gamma_oil=0.8, pb_Mpa=20, co_1Mpa=0.002,
                                    rs_m3m3=300, bo_m3m3=1):
    """
    расчет критического расхода нефтегазовой смеси через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    t_C - температура, С
    d2_mm - диаметр штуцера, мм
    gamma_g - относительная плотность газа по воздуху
    fg - весовая доля газа в потоке, д.ед.
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3

    if fg <= 0: fg = 0.00001
    fo = 1 - fg  # весовая доля нефти в потоке
    m_m = 16  # молекулярный вес пара, моль газа
    z = fl.z  # коэффициент сжимаемости газа
    cvo = 0.24 * 778.169 * PVT.Cvo_kJkgK  # удельная теплоемкость нефти, кДж/кг/К
    cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    roo = roo_kgm3 * 0.062428  # плотность нефти  Ibm/ft3
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    f = cpg / cvg   # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * fo / roo
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fo * cvo) / (fg * cvg + fo * cvo)
    lambd = fg + ((fg * cvg + fo * cvo) * m_m / (z * r_ftIbtlbmmolr))

    def qw(pr2):
        a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
        b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                    fg / n * pr2 ** (-(1 + n) / n))
        c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                    (fg * pr2 ** (-1 / n)) + alf) ** 2
        d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr2 ** (-1 / n)) + alf)
        e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
        return a * (b + c) - d * e

    pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
    pr = pr1
    if pr > 1: pr = 1
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
                (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia/ v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def p_choke_up_gaswater_Perkins_atm(w_kghr, p2_atm, fg, gamma_g=0.55, t_C=20, d2_mm=5):
    """
    расчет давления водогазовой смеси перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    gamma_g - относительная плотность газа по воздуху
    xg - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    """

    def w2(p1_atm):
        return w_choke_gaswater_Perkins_kghr(p1_atm, p2_atm, fg, gamma_g, t_C, d2_mm) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_up_gasoil_Perkins_atm(w_kghr, p2_atm, fg, gamma_g=0.55, gamma_oil=0.8, t_C=20, d2_mm=5, pb_Mpa=20,
                                  co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1):
    """

    расчет давления водогазовой смеси перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    fg - массовая доля газа, д.ед.
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(p1_atm):
        return w_choke_gasoil_Perkins_kghr(p1_atm, p2_atm, fg, gamma_g, t_C, d2_mm, gamma_oil, pb_Mpa,
                                           co_1Mpa, rs_m3m3, bo_m3m3) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_gaswater_Perkins_atm(w_kghr, p1_atm, fg, gamma_g=0.55, t_C=20, d2_mm=5):
    """
    расчет давления водогазовой смеси перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    fg - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    """

    def w2(p2_atm):
        return w_choke_gaswater_Perkins_kghr(p1_atm, p2_atm, fg, gamma_g, t_C, d2_mm) - w_kghr

    p2 = opt.fsolve(w2, p1_atm - 0.00001)
    if w_kghr > w_choke_critgaswater_Perkins_kghr(p1_atm, fg, gamma_g, t_C, d2_mm):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        p2 = 0
    return p2


def p_choke_down_gasoil_Perkins_atm(w_kghr, p1_atm, fg, gamma_g=0.55, gamma_oil=0.8, t_C=20, d2_mm=5, pb_Mpa=20, 
                                    co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет давления газонефтяной смеси перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    fg - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(p2_atm):
        return w_choke_gasoil_Perkins_kghr(p1_atm, p2_atm, fg, gamma_g, t_C, d2_mm, gamma_oil, pb_Mpa,
                                           co_1Mpa, rs_m3m3, bo_m3m3) - w_kghr

    p2 = opt.fsolve(w2, p1_atm - 0.00001)
    if w_kghr > w_choke_critgasoil_Perkins_kghr(p1_atm, fg, gamma_g, t_C, d2_mm, gamma_oil, pb_Mpa,
                                                co_1Mpa, rs_m3m3, bo_m3m3):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        p2 = 0
    return p2


def f_choke_difpressure_gaswater_Perkins_atm(w_kghr, p1_atm, fg, gamma_g=0.55, t_C=20, d2_mm=5):
    """
    расчет разницы давлений до и после штуцера водогазовой смеси по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    fg - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    """

    def w2(f):
        return w_choke_gaswater_Perkins_kghr(p1_atm, p1_atm - f, fg, gamma_g, t_C, d2_mm) - w_kghr

    f = opt.fsolve(w2, 0)
    if w_kghr > w_choke_critgaswater_Perkins_kghr(p1_atm, fg, gamma_g, t_C, d2_mm):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        f = 0
    return f


def f_choke_difpressure_gasoil_Perkins_atm(w_kghr, p1_atm, fg, gamma_g=0.55, gamma_oil=0.8, t_C=20, d2_mm=5, pb_Mpa=20,
                                           co_1Mpa=0.002, rs_m3m3=300, bo_m3m3=1):
    """
    расчет разницы давлений до и после штуцера нефтегазовой смеси по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    fg - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    """

    def w2(f):
        return w_choke_gasoil_Perkins_kghr(p1_atm, p1_atm - f, fg, gamma_g, t_C, d2_mm, gamma_oil, pb_Mpa,
                                           co_1Mpa, rs_m3m3, bo_m3m3) - w_kghr

    f = opt.fsolve(w2, 0)
    if w_kghr > w_choke_critgasoil_Perkins_kghr(p1_atm, fg, gamma_g, t_C, d2_mm, gamma_oil, pb_Mpa,
                                                co_1Mpa, rs_m3m3, bo_m3m3):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        f = 0
    return f


def w_choke_gasoilwater_Perkins_kghr(p1_atm, p2_atm, fg, fo, fw, gamma_g=0.55, t_C=20, d2_mm=5, gamma_oil=0.8, 
                                     pb_Mpa=20, rs_m3m3=300, bo_m3m3=1):
    """
    расчет расхода трехфазной смеси через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    t_C - температура, С
    d2_mm - диаметр штуцера, мм
    gamma_g - относительная плотность газа по воздуху
    fg - весовая доля газа в потоке, д.ед.
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    fo - весовая доля нефти в потоке
    fw -весовая доля воды в потоке
    """

    if fg <= 0: fg = 0.00001
    if fo <= 0: fo = 0.00001
    if fw <= 0: fw = 0.00001
    m_m = 16  # молекулярный вес, моль газа


    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    z = fl.z
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3
    row_kgm3 = steamTable.rhoL_p(p1_atm)  # плотность воды   кг/м3
    cvw = 0.24 * 778.169 * steamTable.CvL_p(p1_atm)  # удельная теплоемкость воды
    # 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
    cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    roo = roo_kgm3 * 0.062428  # плотность нефти  Ibm/ft3
    row = row_kgm3 * 0.062428  # плотность воды  Ibm/ft3
    ro_o = 141.5 / roo_kgm3 / 1000 - 131.5  # плотность нефти, API
    t_F = t_C * 1.8 + 32  # температура, F
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    cvo = 778 * ((0.355 + 0.00176 * ro_o) + (0.0051 + 1.167 * ro_o / 100000) * t_F)  # теплоемкость при постоянном 
    # объеме для нефти (ft-Ibf)/(lbm- OF)
    f = cpg / cvg  # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    p2_psia = p2_atm * 14.2233  # давление за штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * (fo / roo + fw / row)
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fo * cvo + fw * cvw) / (fg * cvg + fo * cvo + fw * cvw)
    lambd = fg + ((fg * cvg + fo * cvo + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
    i = 0
    pr1 = 0.1
    pr11 = 0.2
    """
    Пересчитываю отношение давлений столько раз, пока разница между следующим и предыдущим не будет очень маленькой.
    Отношение давлений служит для того, чтобы найти давление в штуцере.
    """
    while abs(pr1 - pr11) > 0.001 and i < 10:  # пересчитываем удельные теплоемкости и давление в штуцере,
        # пока они не будут соответствовать A*(B+C)=D*E
        def qw(pr2):
            a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
            b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        fg / n * pr2 ** (-(1 + n) / n))
            c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                        (fg * pr2 ** (-1 / n)) + alf) ** 2
            d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        (fg * pr2 ** (-1 / n)) + alf)
            e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
            return a * (b + c) - d * e

        pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
        p3_psia = p1_psia * pr1  # находим давление в штуцере
        pmid_atm = (p1_psia + p3_psia) / 2 / 14.2233  # среднее давление в штуцере и перед штуцером, 
        # пересчитываем удельные теплоемкости для этого давления
        cvw = 0.24 * 778.169 * steamTable.CvL_p(pmid_atm)  # удельная теплоемкость 
        # воды 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
        cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
        cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
        f = cpg / cvg
        n = (fg * f * cvg + fo * cvo + fw * cvw) / (fg * cvg + fo * cvo + fw * cvw)
        lambd = fg + ((fg * cvg + fo * cvo + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
        pr11 = opt.fsolve(qw, 0.01)
        i = i + 1
    if p3_psia > p2_psia:
        pr = p3_psia / p1_psia
    else:
        pr = p2_psia / p1_psia
    if pr > 1: pr = 1
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
                (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia / v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def w_choke_critgasoilwater_Perkins_kghr(p1_atm, fg, fo, fw, gamma_g=0.55, t_C=20, d2_mm=5, gamma_oil=0.8, pb_Mpa=20,
                                         rs_m3m3=300, bo_m3m3=1):
    """
    расчет расхода трехфазной смеси через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    t_C - температура, С
    d2_mm - диаметр штуцера, мм
    gamma_g - относительная плотность газа по воздуху
    fg - весовая доля газа в потоке, д.ед.
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    fo - весовая доля нефти в потоке
    fw -весовая доля воды в потоке
    """

    fl = PVT.FluidMcCain()
    fl.gamma_gas = gamma_g
    fl.gamma_oil = gamma_oil
    fl._bo_m3m3 = bo_m3m3
    fl._pb_bar = pb_Mpa / 10
    fl._rs_m3m3 = rs_m3m3
    fl.calc(p1_atm, t_C)
    z = fl.z
    rog_kgm3 = fl.rho_gas_kgm3
    roo_kgm3 = fl.rho_oil_kgm3


    if fg <= 0: fg = 0.00001
    if fo <= 0: fo = 0.00001
    if fw <= 0: fw = 0.00001
    m_m = 16  # молекулярный вес, моль газа
    row_kgm3 = steamTable.rhoL_p(p1_atm)  # плотность воды   кг/м3
    cvw = 0.24 * 778.169 * steamTable.CvL_p(p1_atm)
    # удельная теплоемкость воды 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
    cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    roo = roo_kgm3 * 0.062428  # плотность нефти  Ibm/ft3
    row = row_kgm3 * 0.062428  # плотность воды  Ibm/ft3
    ro_o = 141.5 / roo_kgm3 / 1000 - 131.5  # плотность нефти, API
    t_F = t_C * 1.8 + 32  # температура, F
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    cvo = 778 * ((0.355 + 0.00176 * ro_o) + (0.0051 + 1.167 * ro_o / 100000) * t_F)  # теплоемкость при постоянном 
    # объеме для нефти (ft-Ibf)/(lbm- OF)
    f = cpg / cvg # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * (fo / roo + fw / row)
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fo * cvo + fw * cvw) / (fg * cvg + fo * cvo + fw * cvw)
    lambd = fg + ((fg * cvg + fo * cvo + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
    i = 0
    pr1 = 0.1
    pr11 = 0.2
    """
    Пересчитываю отношение давлений столько раз, пока разница между следующим и предыдущим не будет очень маленькой.
    Отношение давлений служит для того, чтобы найти давление в штуцере.
    """
    while abs(pr1 - pr11) > 0.001 and i < 10:  # пересчитываем удельные теплоемкости и давление в штуцере, 
        # пока они не будут соответствовать A*(B+C)=D*E
        def qw(pr2):
            a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
            b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        fg / n * pr2 ** (-(1 + n) / n))
            c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                        (fg * pr2 ** (-1 / n)) + alf) ** 2
            d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        (fg * pr2 ** (-1 / n)) + alf)
            e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
            return a * (b + c) - d * e

        pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
        p3_psia = p1_psia * pr1  # находим давление в штуцере
        pmid_atm = (p1_psia + p3_psia) / 2 / 14.2233  # среднее давление в штуцере и перед штуцером,
        #  пересчитываем удельные теплоемкости для этого давления
        cvw = 0.24 * 778.169 * steamTable.CvL_p(
            pmid_atm)  # удельная теплоемкость воды 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
        cvg = 0.24 * 778.169 * PVT.Cvg_kJkgK  # удельная теплоемкость газа при постоянном объеме, кДж/кг/К
        cpg = 0.24 * 778.169 * PVT.Cpg_kJkgK  # удельная теплоемкость газа при постоянном давлении, кДж/кг/К
        f = cpg / cvg
        n = (fg * f * cvg + fo * cvo + fw * cvw) / (fg * cvg + fo * cvo + fw * cvw)
        lambd = fg + ((fg * cvg + fo * cvo + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
        pr11 = opt.fsolve(qw, 0.01)
        i = i + 1
    pr = p3_psia / p1_psia
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
                (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia / v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def p_choke_up_gasoilwater_Perkins_atm(w_kghr, p2_atm, fg, fo, fw, gamma_g=0.55, gamma_oil=0.8, t_C=20, d2_mm=5,
                                       pb_Mpa=20, rs_m3m3=300, bo_m3m3=1):
    """

    расчет давления трехфазной смеси перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    gamma_g - относительная плотность газа по воздуху
    x - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    fg - массовая доля газа, д.ед.
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    fo - весовая доля нефти в потоке
    fw -весовая доля воды в потоке
    """

    def w2(p1_atm):
        return w_choke_gasoilwater_Perkins_kghr(p1_atm, p2_atm, fg, fo, fw, gamma_g, t_C, d2_mm, gamma_oil,
                                                pb_Mpa, rs_m3m3, bo_m3m3) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_gasoilwater_Perkins_atm(w_kghr, p1_atm, fg, fo, fw, gamma_g=0.55, gamma_oil=0.8, t_C=20, d2_mm=5,
                                         pb_Mpa=20, rs_m3m3=300, bo_m3m3=1):
    """
    расчет давления трехфазной смеси перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    fg - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    fo - весовая доля нефти в потоке
    fw -весовая доля воды в потоке
    """

    def w2(p2_atm):
        return w_choke_gasoilwater_Perkins_kghr(p1_atm, p2_atm, fg, fo, fw, gamma_g, t_C, d2_mm, gamma_oil,
                                                pb_Mpa, rs_m3m3, bo_m3m3) - w_kghr

    p2 = opt.fsolve(w2, p1_atm - 0.00001)
    if w_kghr > w_choke_critgasoilwater_Perkins_kghr(p1_atm, fg, fo, fw, gamma_g, t_C, d2_mm, gamma_oil,
                                                     pb_Mpa, rs_m3m3, bo_m3m3):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        p2 = 0
    return p2


def f_choke_difpressure_gasoilwater_Perkins_atm(w_kghr, p1_atm, fg, fo, fw, gamma_g=0.55, gamma_oil=0.8, t_C=20, 
                                                d2_mm=5, pb_Mpa=20, rs_m3m3=300, bo_m3m3=1):
    """
    расчет разницы давлений до и после штуцера трехфазной смеси по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    gamma_g - относительная плотность газа по воздуху
    fg - массовая доля газа, д.ед.
    t_C-  температура ,С
    d2_mm - диаметр штуцера, мм
    gamma_oil - относительная плонтость нефти по воде
    pb_Mpa- давление насыщения, МПа
    co_1Mpa- коэффициент термической сжимаемости нефти 1/MPa
    rs_m3m3 - газонасыщенность нефти м3/м3
    bo_m3m3 - объемный коэффициент  нефти м3/м3
    fo - весовая доля нефти в потоке
    fw -весовая доля воды в потоке
    """

    def w2(f):
        return w_choke_gasoilwater_Perkins_kghr(p1_atm, p1_atm - f, fg, fo, fw, gamma_g, t_C, d2_mm, gamma_oil,
                                                pb_Mpa, rs_m3m3, bo_m3m3) - w_kghr

    f = opt.fsolve(w2, 0)
    if w_kghr > w_choke_critgasoilwater_Perkins_kghr(p1_atm, fg, fo, fw, gamma_g, t_C, d2_mm, gamma_oil,
                                                     pb_Mpa, rs_m3m3, bo_m3m3):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        f = 0
    return f


def w_choke_steam_Chien_kghr(p1_atm, p2_atm, x, d0_mm=5, d1_mm=100):
    """
    расчет расхода насыщенного пара через штуцер по методике Чиена (Миллера)
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    x - cухость пара, д.ед.
    """
    c0 = 0.8  # коэффициент разряда штуцера
    y1 = 0.89487  # коэффициент расширения пара,считают его постоянным,зависит от давлений и диаметров
    a = 0.99998  # эмпирические коэффициенты, рассчитаны из графиков по статье
    b = 1.38
    vfg_m3kg = steamTable.vV_p(p1_atm)  # удельный объем пара при заданном давлении м3/кг
    vf_m3kg = steamTable.vL_p(p1_atm)  # удельный объем жидкости при заданном давлении м3/кг
    vexp_m3kg = a * vfg_m3kg * x ** b + vf_m3kg  # удельный объем смеси м3/кг
    ro_kgm3 = 1 / vexp_m3kg  # плотность смеси кг/м3
    beta = d0_mm / d1_mm
    f_atm = p1_atm - p2_atm
    f_kpa = f_atm * 100
    w_kgsec = 3.512407 * 10 ** (-5) * c0 * y1 * d0_mm ** 2 * (f_kpa * ro_kgm3) ** 0.5 / (1 - beta ** 4) ** 0.5
    w_kghr= w_kgsec * 3600
    return w_kghr


def p_choke_up_steam_Chien_atm(w_kghr, p2_atm, x, d0_mm=5, d1_mm=100):
    """
    расчет давления насыщенного пара перед клапаном по методике Чиена (Миллера)
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    x - cухость пара, д.ед.
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    """

    def w2(p1_atm):
        return w_choke_steam_Chien_kghr(p1_atm, p2_atm, x, d0_mm, d1_mm) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_steam_Chien_atm(w_kghr, p1_atm, x, d0_mm=5, d1_mm=100):
    """
    расчет давления на выходе от давления на входе при разных расходах по методике Чиена (Миллера)
    p1_atm - давление на входе в штуцер, атм
    w_kghr- расход пара, кг/час
    x -  сухость пара, д.ед.
    d0_mm - диаметр штуцера, мм
    d1_mm - диаметр трубы, мм
    """
    c0 = 0.8  # коэффициент разряда штуцера
    y1 = 0.89487  # коэффициент расширения пара
    a = 0.99998  # эмпирические коэффициенты, рассчитаны из графиков по статье
    b = 1.38
    vfg_m3kg = steamTable.vV_p(p1_atm)  # удельный объем пара при заданном давлении м3/кг
    vf_m3kg = steamTable.vL_p(p1_atm)  # удельный объем жидкости при заданном давлении м3/кг
    vexp_m3kg = a * vfg_m3kg * x ** b + vf_m3kg  # удельный объем смеси м3/кг
    ro_kgm3 = 1 / vexp_m3kg  # плотность смеси кг/м3
    beta = d0_mm / d1_mm
    w_kgsec = w_kghr / 3600
    f_kpa = w_kgsec ** 2 * (1 - beta ** 4) * 10 ** 10 / (12.337 * c0 ** 2 * y1 ** 2 * d0_mm ** 4 * ro_kgm3)
    p2_atm = p1_atm - f_kpa / 100
    return p2_atm


def f_choke_difpressure_steam_Chien_atm(w_kghr, p1_atm, x, d0_mm=5, d1_mm=100):
    """
    расчет разницы давлений до и после штуцера насыщенного пара по методике Чиена (Миллера)
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    x - cухость пара, д.ед.
    d0_mm - диаметр штуцера, мм
    d1_mm -  диаметр трубы, мм
    """

    def w2(f):
        return w_choke_steam_Chien_kghr(p1_atm, p1_atm - f, x, d0_mm, d1_mm) - w_kghr

    return opt.fsolve(w2, 0)


def w_choke_steam_AlSafran_kghr(p1_atm, p2_atm, xg, dchoke_mm=5):
    """
    расчет расхода насыщенного пара через штуцер по методике Альсафран
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    dchoke_mm - диаметр штуцера,мм
    xg - массовая доля пара, д.ед.
    """
    if xg <= 0: xg = 0.0001
    cd = 0.75  # дисчарч коэффициент для штуцера (коэффициент разряда штуцера)
    g_c = 9.8  # ускорение свободного падения
    p1_pa = p1_atm * 100  # давление до штуцера, Па
    p2_pa = p2_atm * 100  # давление после штуцера, Па
    vgl = steamTable.vV_p(p1_atm)  # удельный объем пара, м3/кг
    vl = steamTable.vL_p(p1_atm)  # удельный объем воды, м3/кг
    cl = steamTable.CvL_p(p1_atm)  # удельная теплоемкость воды, кДж/кг/К
    cvg = steamTable.CvV_p(p1_atm)  # удельная теплоемкость пара при постоянном объеме, кДж/кг/К
    cpg = steamTable.CpV_p(p1_atm)  # удельная теплоемкость пара при постоянном давлении, кДж/кг/К
    k = cpg / cvg # показатель адиабаты
    d2 = dchoke_mm / 1000  # диаметр клапана, м
    r = (1 + xg * (vgl / vl - 1)) ** 0.5 * (1 + 0.6 / 2.718 ** (5 * xg))
    alfa = r * (1 - xg) * vl / (xg * vgl)
    a2 = 3.14 * d2 ** 2 / 4
    n = (xg * k * cvg + (1 - xg) * cl) / (xg * cvg + (1 - xg) * cl)

    def y(x):
        y = (alfa * (1 - x) + n / (n - 1)) / (n / (n - 1) + n / 2 * (1 + alfa * x ** (1 / n)) ** 2) - x ** (1 - 1 / n)
        return (y)

    rc = opt.fsolve(y, 0.5)
    r= p2_pa / p1_pa
    if rc >= r:
        pr = rc
    else:
        pr = r
    if pr > 1: pr = 1
    ab = n / (n - 1) * (1 - pr ** ((n - 1) / n)) + alfa * (1 - pr)
    ac = xg * vgl * (pr ** (-1 / n) + alfa) ** 2 * (xg + 1 / r * (1 - xg))
    wi = a2 * (288 * g_c * cd ** 2 * p1_pa * ab / ac) ** 0.5 * 3600
    return wi


def w_choke_critsteam_AlSafran_kghr(p1_atm, xg, dchoke_mm=5):
    """
    расчет критического расхода насыщенного пара через штуцер по методике Альсафран
    p1_atm - давление на входе в штуцер, атм
    dchoke_mm - диаметр штуцера,мм
    xg - массовая доля пара, д.ед.
    """
    if xg <= 0: xg = 0.0001
    cd = 0.75  # дисчарч коэффициент для штуцера (коэффициент разряда штуцера)
    g_c = 9.8  # ускорение свободного падения
    p1_pa = p1_atm * 100  # давление до штуцера, Па
    vgl = steamTable.vV_p(p1_atm)  # удельный объем пара, м3/кг
    vl = steamTable.vL_p(p1_atm)  # удельный объем воды, м3/кг
    cl = steamTable.CvL_p(p1_atm)  # удельная теплоемкость воды, кДж/кг/К
    cvg = steamTable.CvV_p(p1_atm)  # удельная теплоемкость пара при постоянном объеме, кДж/кг/К
    cpg = steamTable.CpV_p(p1_atm)  # удельная теплоемкость пара при постоянном давлении, кДж/кг/К
    k = cpg / cvg   # показатель адиабаты
    d2 = dchoke_mm / 1000  # диаметр клапана, м
    r = (1 + xg * (vgl / vl - 1)) ** 0.5 * (1 + 0.6 / 2.718 ** (5 * xg))
    alfa = r * (1 - xg) * vl / (xg * vgl)
    a2 = 3.14 * d2 ** 2 / 4
    n = (xg * k * cvg + (1 - xg) * cl) / (xg * cvg + (1 - xg) * cl)

    def y(x):
        y = (alfa * (1 - x) + n / (n - 1)) / (n / (n - 1) + n / 2 * (1 + alfa * x ** (1 / n)) ** 2) - x ** (1 - 1 / n)
        return (y)

    rc = opt.fsolve(y, 0.5)
    pr = rc
    if pr > 1: pr = 1
    ab = n / (n - 1) * (1 - pr ** ((n - 1) / n)) + alfa * (1 - pr)
    ac = xg * vgl * (pr ** (-1 / n) + alfa) ** 2 * (xg + 1 / r * (1 - xg))
    wi = a2 * (288 * g_c * cd ** 2 * p1_pa * ab / ac) ** 0.5 * 3600
    return wi


def p_choke_up_steam_AlSafran_atm(w_kghr, p2_atm, xg, dchoke_mm=5):
    """
    расчет давления насыщенного пара перед клапаном по методике Альсафран
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    xg - cухость пара, д.ед.
    dchoke_mm - диаметр штуцера, мм
    """

    def w2(p1_atm):
        return w_choke_steam_AlSafran_kghr(p1_atm, p2_atm, xg, dchoke_mm) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_steam_AlSafran_atm(w_kghr, p1_atm, xg, dchoke_mm=5):
    """
    расчет давления насыщенного пара перед клапаном по методике Альсафран
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    xg - cухость пара, д.ед.
    dchoke_mm - диаметр штуцера, мм
    """

    def w2(p2_atm):
        return w_choke_steam_AlSafran_kghr(p1_atm, p2_atm, xg, dchoke_mm) - w_kghr

    p2 = opt.fsolve(w2, p1_atm - 0.00001)
    if w_kghr > w_choke_critsteam_AlSafran_kghr(p1_atm, xg):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        p2 = 0
    return p2


def f_choke_difpressure_steam_AlSafran_atm(w_kghr, p1_atm, xg, dchoke_mm=5):
    """
    расчет разницы давлений до и после штуцера насыщенного пара по методике Aльсафран
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    xg - cухость пара, д.ед.
    dchoke_mm - диаметр штуцера, мм
    """

    def w2(f):
        return w_choke_steam_AlSafran_kghr(p1_atm, p1_atm - f, xg, dchoke_mm) - w_kghr

    f = opt.fsolve(w2, 0)
    if w_kghr > w_choke_critsteam_AlSafran_kghr(p1_atm, xg):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        f = 0
    return f


def w_choke_steam_Perkins_kghr(p1_atm, p2_atm, fg, d2_mm=5):
    """
    расчет расхода насыщенного пара через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    t_C - температура, С
    d2_mm - диаметр штуцера, мм
    fg - весовая доля газа в потоке или по другому сухость пара, д.ед.
    """
    if fg <= 0: fg = 0.00001
    fw = 1 - fg  # весовая доля воды в потоке
    m_m = 18  # молекулярный вес пара, моль пара
    z = 0.999  # коэффициент сжимаемости пара
    rog_kgm3 = steamTable.rhoV_p(p1_atm)  # плотность пара кг/м3
    row_kgm3 = steamTable.rhoL_p(p1_atm)  # плотность воды   кг/м3
    cvw = 0.24 * 778.169 * steamTable.CvL_p(
        p1_atm)  # удельная теплоемкость воды 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
    cvg = 0.24 * 778.169 * steamTable.CvV_p(p1_atm)  # удельная теплоемкость пара при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * steamTable.CpV_p(p1_atm)  # удельная теплоемкость пара при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    row = row_kgm3 * 0.062428  # плотность воды  Ibm/ft3
    # t_F=t_C*1.8 + 32 # температура, F
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    f = cpg / cvg # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    p2_psia = p2_atm * 14.2233  # давление за штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * fw / row
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fw * cvw) / (fg * cvg+ fw * cvw)
    lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
    i = 0
    pr1 = 0.1
    pr11 = 0.2
    """
    Пересчитываю отношение давлений столько раз, пока разница между следующим и предыдущим не будет очень маленькой.
    Отношение давлений служит для того, чтобы найти давление в штуцере.
    """
    while abs(pr1 - pr11) > 0.001 and i < 10:  # пересчитываем удельные теплоемкости и давление в штуцере,
        #  пока они не будут соответствовать A*(B+C)=D*E
        def qw(pr2):
            a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
            b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        fg / n * pr2 ** (-(1 + n) / n))
            c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                        (fg * pr2 ** (-1 / n)) + alf) ** 2
            d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        (fg * pr2 ** (-1 / n)) + alf)
            e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
            return a * (b + c) - d * e

        pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
        p3_psia = p1_psia * pr1  # находим давление в штуцере
        pmid_atm = (p1_psia + p3_psia) / 2 / 14.2233  # среднее давление в штуцере и перед штуцером,
        #  пересчитываем удельные теплоемкости для этого давления
        cvw = 0.24 * 778.169 * steamTable.CvL_p(pmid_atm)  # удельная теплоемкость жидкости
        #  1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
        cvg = 0.24 * 778.169 * steamTable.CvV_p(pmid_atm)  # удельная теплоемкость пара при постоянном объеме, кДж/кг/К
        cpg = 0.24 * 778.169 * steamTable.CpV_p(
            pmid_atm)  # удельная теплоемкость пара при постоянном давлении, кДж/кг/К
        f = cpg / cvg
        n = (fg * f * cvg + fw * cvw) / (fg * cvg + fw * cvw)
        lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
        pr11 = opt.fsolve(qw, 0.01)
        i = i + 1
    if p3_psia > p2_psia:
        pr = p3_psia / p1_psia
    else:
        pr = p2_psia / p1_psia
    if pr > 1: pr = 1
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
                (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia / v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def w_choke_critsteam_Perkins_kghr(p1_atm, fg, d2_mm=5):
    """
    расчет критического расхода насыщенного пара через штуцер по методике Перкинса
    p1_atm - давление на входе в штуцер, атм
    p2_atm - давление на выходе из штуцера, атм
    t_C - температура, С
    d2_mm - диаметр штуцера, мм
    fg - весовая доля газа в потоке или по другому сухость пара, д.ед.
    """
    if fg <= 0: fg = 0.00001
    fw = 1 - fg  # весовая доля воды в потоке
    m_m = 18  # молекулярный вес пара, моль пара
    z = 0.999  # коэффициент сжимаемости пара
    rog_kgm3 = steamTable.rhoV_p(p1_atm)  # плотность пара кг/м3
    row_kgm3 = steamTable.rhoL_p(p1_atm)  # плотность воды   кг/м3
    cvw = 0.24 * 778.169 * steamTable.CvL_p(
        p1_atm)  # удельная теплоемкость воды 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
    cvg = 0.24 * 778.169 * steamTable.CvV_p(p1_atm)  # удельная теплоемкость пара при постоянном объеме, кДж/кг/К
    cpg = 0.24 * 778.169 * steamTable.CpV_p(p1_atm)  # удельная теплоемкость пара при постоянном давлении, кДж/кг/К
    d1_mm = 100  # диаметр трубы до штуцера мм
    g_c = 32.17  # (lbm-ft)/(lbf-second^2)
    row = row_kgm3 * 0.062428  # плотность воды  Ibm/ft3
    r_ftIbtlbmmolr = 1545.348  # универсальная газовая постоянная  (ft-Ibf)/(lbm mol-R)
    f = cpg / cvg # показатель адиабаты
    p1_psia = p1_atm * 14.2233  # давление перед штуцером в psia
    d1_ft = d1_mm / 304.8  # диаметр трубы в ft
    d2_ft = d2_mm / 304.8  # диаметр штуцера в ft
    v1_ft3lbm = 16.01845 / rog_kgm3  # удельный объем газа ft3/Ibm
    alf = (1 / v1_ft3lbm) * fw / row
    a1_ft2 = 3.14 * d1_ft ** 2 / 4
    a2_ft2 = 3.14 * d2_ft ** 2 / 4
    n = (fg * f * cvg + fw * cvw) / (fg * cvg + fw * cvw)
    lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
    i = 0
    pr1 = 0.1
    pr11 = 0.2
    """
    Пересчитываю отношение давлений столько раз, пока разница между следующим и предыдущим не будет очень маленькой.
    Отношение давлений служит для того, чтобы найти давление в штуцере.
    """
    while abs(
            pr1 - pr11) > 0.001 and i < 10:  # пересчитываем удельные теплоемкости и давление в штуцере,
        # пока они не будут соответствовать A*(B+C)=D*E
        def qw(pr2):
            a = 2 * lambd * (1 - pr2 ** ((n - 1) / n)) + 2 * alf * (1 - pr2)
            b = (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        fg / n * pr2 ** (-(1 + n) / n))
            c = (a2_ft2 / a1_ft2) ** 2 * fg / n * (fg + alf) ** 2 * pr2 ** (-(1 + n) / n) / (
                        (fg * pr2 ** (-1 / n)) + alf) ** 2
            d = (1 - (a2_ft2 / a1_ft2) * ((fg + alf) / ((fg * pr2 ** (-1 / n)) + alf)) ** 2) * (
                        (fg * pr2 ** (-1 / n)) + alf)
            e = lambd * (n - 1) / n * pr2 ** (-1 / n) + alf
            return a * (b + c) - d * e

        pr1 = opt.fsolve(qw, 0.01)  # находим отношение давления в штуцере и давления перед штуцером
        p3_psia = p1_psia * pr1  # находим давление в штуцере
        pmid_atm = (p1_psia + p3_psia) / 2 / 14.2233  # среднее давление в штуцере и перед штуцером,
        #  пересчитываем удельные теплоемкости для этого давления
        cvw = 0.24 * 778.169 * steamTable.CvL_p(pmid_atm)  # удельная теплоемкость жидкости
        # 1 кДж/кгК= 1*0,24*778.169=(ft-Ibf)/(lbm- OF)
        cvg = 0.24 * 778.169 * steamTable.CvV_p(pmid_atm)  # удельная теплоемкость пара при постоянном объеме, кДж/кг/К
        cpg = 0.24 * 778.169 * steamTable.CpV_p(pmid_atm)  # удельная теплоемкость пара при постоянном давлении,кДж/кг/К
        f = cpg / cvg
        n = (fg * f * cvg + fw * cvw) / (fg * cvg + fw * cvw)
        lambd = fg + ((fg * cvg + fw * cvw) * m_m / (z * r_ftIbtlbmmolr))
        pr11 = opt.fsolve(qw, 0.01)
        i = i + 1
    pr = p3_psia / p1_psia
    if pr > 1: pr = 1
    ab = (lambd * (1 - pr ** ((n - 1) / n)) + alf * (1 - pr)) / (
                (1 - ((a2_ft2 / a1_ft2) ** 2) * ((fg + alf) / ((fg * pr ** (-1 / n)) + alf)) ** 2) * (
                    (fg * pr ** (-1 / n)) + alf) ** 2)
    wi_lbmsec = a2_ft2 * ((288 * g_c * p1_psia / v1_ft3lbm) * ab) ** 0.5
    wi_kghr = wi_lbmsec * 0.45359 * 3600  # перевод в кг/час
    return wi_kghr


def p_choke_up_steam_Perkins_atm(w_kghr, p2_atm, fg, d2_mm=5):
    """
    расчет давления насыщенного пара перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p2_atm - давление на выходе из штуцера, атм
    fg - cухость пара, д.ед.
    d2_mm - диаметр штуцера, мм
    """

    def w2(p1_atm):
        return w_choke_steam_Perkins_kghr(p1_atm, p2_atm, fg, d2_mm) - w_kghr

    return opt.fsolve(w2, p2_atm)


def p_choke_down_steam_Perkins_atm(w_kghr, p1_atm, fg, d2_mm=5):
    """
    расчет давления насыщенного пара перед клапаном по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    fg - cухость пара, д.ед.
    d2_mm - диаметр штуцера, мм
    """

    def w2(p2_atm):
        return w_choke_steam_Perkins_kghr(p1_atm, p2_atm, fg, d2_mm) - w_kghr

    p2 = opt.fsolve(w2, p1_atm - 0.00001)
    if w_kghr > w_choke_critsteam_Perkins_kghr(p1_atm, fg):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        p2 = 0
    return p2


def f_choke_difpressure_steam_Perkins_atm(w_kghr, p1_atm, fg, d2_mm=5):
    """
    расчет разницы давлений до и после штуцера насыщенного пара по методике Перкинса
    w_kghr- расход насыщенного пара, кг/час
    p1_atm - давление на входе в штуцер, атм
    fg - cухость пара, д.ед.
    d2_mm - диаметр штуцера, мм
    """

    def w2(f):
        return w_choke_steam_Perkins_kghr(p1_atm, p1_atm - f, fg, d2_mm) - w_kghr

    f = opt.fsolve(w2, 0)
    if w_kghr > w_choke_critsteam_Perkins_kghr(p1_atm, fg):
        print('Ошибка! Слишком большой расход для такого входного давления!')
        f = 0
    return f

