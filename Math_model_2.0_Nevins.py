
import math
import random
import time

with open('data.txt', 'r') as f:
    data_ = eval(f.read())

with open('T_ion.txt', 'r') as g:
    T_ion_list = eval(g.read())

#-------------------------Константы--------------------------
energy_He1 = 4.5                    #Энергия первой альфа частицы МэВ
energy_He2 =2.09                    #Энергия второй и третей альфа частицы МэВ
weight_electron = 9.109e-31         #Масса элетрона [кг]
weight_B = 1.82e-26                 #Масса бора  [кг]
weight_He = 6.64e-27                #Масса Гелия 
weight_H = 1.6e-27                  #Масса Водорода (протона)
k = 1.38e-23                        #Постоянная Больцмана [Дж]
charge_H = 1                        #Заряд Водорода (протона)
charge_B = 5                        #Заряд Бора 
charge_He = 2                       #Заряд гелия 
energy_fusion = 1.389e-12           #Энергия слияния
e_0 = 8.85e-12                      #электрическая постоянная 
charge_electron = 1.6e-19           #Заряд электрона 
log_K = 16                          #Кулоновский логарифм
C_f = 8                             #Какая то имерическая константа 
h_ion = 1                           #Доля мощности внешнего нагрева для иона
h_electron = 0                      #Доля мощности внешнего нагрева для электронов

#-------Имерические коэффициенты для скорости реакции-------
P_1 = 4.4467e-14              
P_2 = -5.9351e-2
P_3 = 2.0165e-1
P_4 = 1.0404e-3
P_5 = 2.7621e-3
P_6 = -9.1653e-6
P_7 = 9.8305e-7
#------------------------------------------------------------

MrC2 = 859526                       #энергия по Эйнштейну, переведенная в кэВ
Energy_G = 22580                    #Энергия Гамова
alpha_struct = 7.297352e-3          #Постоянная тонкой структуры
radius_electron = 2.8179e-15        #Радиус электрона 
speed_light = 299792458             #Скорость света в вакууме

C_b = alpha_struct*(radius_electron**2)*weight_electron*(speed_light**3)        #Какая то постоянная для тормозного излучения
#------------------------------------------------------------


def prepar_calculate (data_val:list):
    '''
        Находит промежуточные вычисления - концентрации
    '''
    concentr_B = data_val[0]*data_val[1]                             #концентрация бора
    concentr_ion = concentr_B+data_val[0]                            #концентрация ионов плазмы
    return (concentr_B,concentr_ion)

def speed_reaction (data_val:list,T_ion:int,concentr_B:float):
    '''
        Считает скорость реакции для данного набора 
                        данных
    '''
    teta = T_ion/(1-T_ion*
         ((P_2+T_ion*(P_4+T_ion*P_6))/(1+T_ion*(P_3+T_ion*(P_5+T_ion*P_7))))     #коэффициент - эквивалент температуре
    )
    ksi = (Energy_G/(4*teta))**(1/3)                                             #промежуточный коэффициент
    speed_react_NR = P_1*teta*((ksi/(MrC2*T_ion**3))**(1/2))*math.exp(-3*ksi)    #Скорость реакции (нерелятивстский случай)
    speed_react_R = 5.41e-21*((1/T_ion)**(3/2))*math.exp (-148/T_ion)            #Скорость реакции (релятивистский случай)
    speed_react = speed_react_NR + speed_react_R                                 #Скорость реакции суммарная
    power_fus = data_val[0]*concentr_B*speed_react*energy_fusion                 #Тепловая мощность реакции
    return (power_fus,speed_react)

def intermed_calculate (data_val:list,concentr_B:float,speed_react: float,):
    '''
        Считает прмежуточные значения, необходимые для 
                     дальнейших вычислений 
    '''
    concentr_He = 3*data_val[0]*concentr_B*speed_react*data_val[2]                                #Концентрация гелия 
    fraction_He = concentr_He/data_val[0]
    concentr_electron = data_val[0]+5*concentr_B+2*concentr_He                                    #Концентрация электронов в плазме
    charge_effect = ((data_val[0]*charge_H**2+concentr_B*charge_B**2+concentr_He*charge_He**2)/   #Квадрат эффективного заряда
                    (data_val[0]*charge_H+concentr_B*charge_B+concentr_He*charge_He))
    return (concentr_electron,charge_effect,fraction_He)
def power_ion_electron (data_val:list,T_ion:int,T_electron:float,concentr_electron:float,concentr_B:float):
    '''
        Считает мощность, передаваемую от ионов 
                     к электронам 
    '''
    time_H_electron = (((3*math.pi*math.sqrt(2*math.pi)*(e_0**2)*weight_H*weight_electron)/
                       ((charge_H**2)*(charge_electron**4)*concentr_electron*log_K))*               #время столкновения протонов с электронами
                       ((k*T_electron)/(weight_electron*8.61e-8))**(3/2)
                      )
    time_B_electron = (((3*math.pi*math.sqrt(2*math.pi)*(e_0**2)*weight_B*weight_electron)/
                       ((charge_B**2)*(charge_electron**4)*concentr_electron*log_K))*               #время столкновения бора с электронами
                       ((k*T_electron)/(weight_electron*8.61e-8))**(3/2)
                      )
    power_ion_electron_ = (((1.5*data_val[0]*k*((T_ion-T_electron)/8.61e-8))/time_H_electron)+
                           ((1.5*concentr_B*k*((T_ion-T_electron)/8.61e-8))/time_B_electron)       #Можность, передаваемая от ионов к электронам
                          )
    return (power_ion_electron_)

def approximation_func (y:float):
    '''
        Аппрокисмирующая функция, нужная для расчета
                 долей энергии заряженных 
        продуктов, передаваемой ионам  и электронам
    '''
    try:
        f = ((4/9)*math.pi*math.sqrt(3)+(1/3)*math.log(abs(1-math.sqrt(y)))-(1/3)*math.log(1+math.sqrt(y))+
         (1/6)*math.log(1-math.sqrt(y)+y)+(1/3)*math.sqrt(3)*math.atan((math.sqrt(3*y))/(-2+math.sqrt(y)))-
         (1/6)*math.log(1+math.sqrt(y)+y)-(1/3)*math.sqrt(3)*math.atan((math.sqrt(3*y))/(2+math.sqrt(y)))-
         (1/3)*math.log(abs(1-y))+(1/6)*math.log(1+y+y**2)-(1/3)*math.sqrt(3)*math.atan((y*math.sqrt(3))/(2+y))
        )   
    except Exception as e:
        y = y +0.01
        f = ((4/9)*math.pi*math.sqrt(3)+(1/3)*math.log(abs(1-math.sqrt(y)))-(1/3)*math.log(1+math.sqrt(y))+
         (1/6)*math.log(1-math.sqrt(y)+y)+(1/3)*math.sqrt(3)*math.atan((math.sqrt(3*y))/(-2+math.sqrt(y)))-
         (1/6)*math.log(1+math.sqrt(y)+y)-(1/3)*math.sqrt(3)*math.atan((math.sqrt(3*y))/(2+math.sqrt(y)))-
         (1/3)*math.log(abs(1-y))+(1/6)*math.log(1+y+y**2)-(1/3)*math.sqrt(3)*math.atan((y*math.sqrt(3))/(2+y))
        )
    return (f)

def alpha_ion_electron (data_val:list,T_electron:float,concentr_electron:float,concentr_B:float):
    '''
        Считает доли энергии заряженных продуктов, передаваемой
                    ионам и электронам 

    '''
    speed_brake = (((0.75*math.sqrt(math.pi)*(weight_electron/concentr_electron)*
                  (((data_val[0]*charge_H**2)/weight_H)+((concentr_B*charge_B**2)/weight_B)))**(1/3))*       #Скорость торможения на ионах
                  ((2*k*T_electron)/(weight_electron*8.61e-8))**0.5
                  )
    energy_brake = 0.5*weight_He*speed_brake**2*6.242e12                     #Энергия торможения на ионах
    y_He1 = energy_brake/energy_He1                                          #Отношение Энергии торможения на ионах к Энергии гелия 1
    y_He2 = energy_brake/energy_He2                                          #Отношение Энергии торможения на ионах к Энергии гелия 2
    delta_energy_He1 = y_He1*approximation_func(y_He1)                       #Промежутоное значение для гелия с энергией 1 
    delta_energy_He2 = y_He2*approximation_func(y_He2)                       #Промежутоное значение для гелия с энергией 2
    if delta_energy_He1 >=1:
        delta_energy_He1 = 1
    if delta_energy_He2 >=1:
        delta_energy_He2 = 1
    alpha_ion = (delta_energy_He1*energy_He1+2*delta_energy_He2*energy_He2)/(energy_He1+2*energy_He2)   #Доля энергии, передаваемой ионам
    alpha_electron = 1-alpha_ion                                                                        #Доля энергии, передаваемой электронам 
    return (alpha_ion,alpha_electron)

def bremsstrahlung (T_electron:float,charge_effect:float,concentr_electron:float):
    '''
        Находит тормозое излучение в плазме
    '''
    temp_rest = (k*T_electron)/(weight_electron*(speed_light**2)*8.61e-8)                          #Температура враженная в еденицах покоя
    power_Br_e_ion = ((32/3)*math.sqrt(2/math.pi)*C_b*(concentr_electron**2)*charge_effect*
                       math.sqrt(temp_rest)*(0.68+0.32*math.exp(-4.4*temp_rest)+2.070*temp_rest)   #Мощность тормозного излучения на ионах
                     )
    power_Br_e_e = (4*C_f*(math.pi)**(-1/2)*C_b*(concentr_electron**2)*(temp_rest**1.5)*
                     (1+0.64*temp_rest+6.6*temp_rest**2-22.6*temp_rest**3+33.8*temp_rest**4-24.7*temp_rest**5+7.1*temp_rest**6) #Мощность ТИ на эл.
                   )
    power_bremsstrahlung = power_Br_e_e+power_Br_e_ion                                              #Суммарная мощность ТИ
    return (power_bremsstrahlung)

def power_extern (data_val:list,concentr_ion:float,T_ion:int,concentr_electron:float,T_electron:float,power_bremsstrahlung:float,power_fus:float):
    '''
        Находит мощность внешнего нагрева плазмы
    '''
    power_ext = (((3/2)*k)/(data_val[2]*8.61e-8))*(concentr_ion*T_ion+concentr_electron*T_electron)+power_bremsstrahlung-power_fus #Мощность вн. нагрева
    return (power_ext)

def equation_0 (data_val:list,alpha_ion:float,power_fus:float,power_ext:float,concentr_ion:float,T_ion,power_ion_electron_:float):
    '''
        Считает конечное уравнение для нахождения 
                    погрешности
    '''
    equation = (alpha_ion*power_fus+h_ion*power_ext-((3/2*concentr_ion*k*(T_ion/(8.61e-8)))/data_val[2])-power_ion_electron_)/1e6
    return (equation)

def T_electron_detect (data_state:list,T_ion:int,T_electron:float):
    concentr_B,concentr_ion = prepar_calculate(data_state)
    power_fus, speed_react = speed_reaction(data_state, T_ion, concentr_B)
    concentr_electron, charge_effect,fraction_He = intermed_calculate (data_state,concentr_B,speed_react)
    power_ion_electron_ = power_ion_electron (data_state,T_ion,T_electron,concentr_electron,concentr_B)
    alpha_ion, alpha_alectron = alpha_ion_electron (data_state,T_electron,concentr_electron,concentr_B)
    power_bremsstrahlung = bremsstrahlung(T_electron,charge_effect,concentr_electron)
    power_ext = power_extern (data_state,concentr_ion,T_ion,concentr_electron,T_electron,power_bremsstrahlung,power_fus)
    equation = equation_0 (data_state,alpha_ion,power_fus,power_ext,concentr_ion,T_ion,power_ion_electron_)
    return equation


def T_electron_main (data_val:list):
    '''
        Итерационно находит температуру электронов 
            с помощью метода, придуманного 
                        мной
    '''
    data_new = []                                                       #Список с состояниями и соотеветствующими им температуры
    for data_state in data_val:
        array_temp = []                                                 #Список температур состояния 
        for T_ion in T_ion_list:
            a = 0                                                       #Левая часть интервала моего метода
            b = 450                                                     #Правая часть интервала моего метода
            T_electron = (b+a)/2                                        #Температура электронов, взятая из моего интервала
            error = T_electron_detect (data_state,T_ion,T_electron)     #Отклонение от 0 (верного значения)
            i = 0                                                       #Просто счетчик итераций
            while ((abs(error))>0.05 and i<20):
                if error > 0:
                    b = T_electron
                    T_electron = (b+a)/2
                else:
                    a = T_electron
                    T_electron = (b+a)/2
                error = T_electron_detect (data_state,T_ion,T_electron) 
                i= i+1   
            array_temp_smal = [T_ion,T_electron]                        #Список температур электронов и ионв (2 значения)
            array_temp.append (array_temp_smal)
        data_new.append([data_state,array_temp])
    return (data_new)
    
def coeff_gain_power_Pfus_Pb_func (data_state:list,T_ion:int, T_electron:float):
    '''
        Находит коэффициент мощности 
           для заданного состояния
    '''
    concentr_B,concentr_ion = prepar_calculate(data_state)
    power_fus, speed_react = speed_reaction(data_state, T_ion, concentr_B)
    concentr_electron, charge_effect, fraction_He = intermed_calculate (data_state,concentr_B,speed_react)
    power_ion_electron_ = power_ion_electron (data_state,T_ion,T_electron,concentr_electron,concentr_B)
    alpha_ion, alpha_alectron = alpha_ion_electron (data_state,T_electron,concentr_electron,concentr_B)
    power_bremsstrahlung = bremsstrahlung(T_electron,charge_effect,concentr_electron)
    coeff_gain_power = power_fus/((((3/2)*k)/(data_state[2]*8.61e-8))*(concentr_ion*T_ion+concentr_electron*T_electron)+power_bremsstrahlung-power_fus)
    P_fus_Pb = power_fus/power_bremsstrahlung
    return (coeff_gain_power,P_fus_Pb,fraction_He)

def coeff_gain_power_main (data_val:list,):
    '''
        Находит для каждого микросостояния 
          коэффициент усиления мощности
            и записывает все в список
    '''
    data_new = []
    for data_state in data_val:
        array_temp_Q_full = []
        for data_temp in data_state[1]:
            coeff_gain_power, Pfus_Pb, fraction_He = coeff_gain_power_Pfus_Pb_func(data_state[0],data_temp[0],data_temp[1])
            array_temp_Q = [data_temp[0],data_temp[1],coeff_gain_power,Pfus_Pb, fraction_He]
            array_temp_Q_full.append(array_temp_Q)
        data_new.append([data_state[0], array_temp_Q_full])
    return (data_new)

if __name__ == '__main__':
    
    data_new = T_electron_main (data_)
                                                                            #Модуль создания списка с готовыми температурами    
    with open('Temp_electron_Nevins.txt', 'w') as f:
        f.write(str(data_new))

    with open('Temp_electron_Nevins.txt', 'r') as f:
        data_ = eval(f.read())
    
    data_new = coeff_gain_power_main (data_)                                #Модуль создания списка с коэффициентом усиления мощности

    with open('T_el__coeff_power__Pfus_pb__XHe_Nevins.txt', 'w') as g:
        g.write(str(data_new))
    print ('end')
    
    
    
