from settings import *


#--------Константы, которые задаем--------
concentr_prot_min = CONCENTR_PROTONS_START      #минимальная Концентрация протонов [1/м^3]
concentr_prot_max = CONCENTR_PROTONS_END        #максимальная Концентрация протонов [1/м^3]
count_concentr_prot = STEP_CONCENTR_PROTONS     #Количество концентраций протонов в заданном диапазоне
share_bor_min = SHARE_BOR_START                 #минимальное Отношение концентрации бора к концентрации протонов (доля бора) 
share_bor_max = SHARE_BOR_END                   #максимальное Отношение концентрации бора к концентрации протонов (доля бора)
count_bor = STEP_SHARE_BOR                      #Количество доль бора в заданном диапозоне 
time_hold_min = TIME_HOLD_START                 #минимальное Время удержания плазмы [с]
time_hold_max = TIME_HOLD_END                   #максимальное Время удержания плазмы [с]
count_time_hold = STEP_TIME_HOLD                #количесвто времен удержания в заданном интервале
approx_Te = APPROXE_T_ELECTRON                  #Начальное приближение температуры электронов [кэВ]
count_T = 46                                    #Количество разбиений по температуре 
#--------------------------------------------

def Creat_DATA():
    data_ = []
    for i in range (count_concentr_prot+1):
        for j in range (count_bor+1):
            for k in range(count_time_hold+1):
                arr_1 = [concentr_prot_min+i*((concentr_prot_max-concentr_prot_min)/count_concentr_prot),
                         share_bor_min+j*((share_bor_max-share_bor_min)/count_bor),
                         time_hold_min+k*((time_hold_max-time_hold_min)/count_time_hold)
                        ]
                data_.append (arr_1)
    return (data_)


T_i = T_ION_START
arr = []
while T_i < T_ION_END+1:
    arr.append(T_i)
    T_i = T_i + STEP_T_ION

data_ = Creat_DATA ()


with open('T_ion.txt', 'w') as f:
    f.write(str(arr))

with open('data.txt', 'w') as f:
    f.write(str(data_))