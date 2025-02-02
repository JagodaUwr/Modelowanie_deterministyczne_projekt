from kod import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as am

alpha = 0.025
ro = 1.25
c = 1005

# dyskretyzacja przestrzeni
hx = 0.5
xmax = 18
ymax = 12
x = np.arange(0, xmax+hx, hx)
y = np.arange(0, ymax+hx, hx)
X, Y = np.meshgrid(x, y)

# dyskretyzacja czasu
Tmax = 12*60*60
ht = 0.5
t = np.arange(0, Tmax, ht)


def find_wall_indices(x, y, left, right, down, up):
    left_wall = np.ravel_multi_index((np.arange(down, up + 1), np.full(up - down + 1, left)), (len(y), len(x)))
    right_wall = np.ravel_multi_index((np.arange(down, up + 1), np.full(up - down + 1, right)), (len(y), len(x)))
    bottom_wall = np.ravel_multi_index((np.full(right - left + 1, down), np.arange(left, right + 1)), (len(y), len(x)))
    top_wall = np.ravel_multi_index((np.full(right - left + 1, up), np.arange(left, right + 1)), (len(y), len(x)))

    return left_wall, right_wall, bottom_wall, top_wall


def find_index(arr, value):
    idx = np.where(arr == value)[0]
    if len(idx) > 0:
        return int(idx[0])
    else:
        raise ValueError(f"Nie znaleziono wartości {value} w tabeli.")


# room 1
left1 = find_index(x, 0.) + 1
right1 = find_index(x, 10.)
down1 = find_index(y, 0.) + 1
up1 = find_index(y, 6.)

left_wall1, right_wall1, bottom_wall1, top_wall1 = find_wall_indices(x, y, left1, right1, down1, up1)

# room 2
left2 = find_index(x, 0.) + 1
right2 = find_index(x, 10.)
down2 = find_index(y, 6.) + 1
up2 =  find_index(y, ymax)

left_wall2, right_wall2, bottom_wall2, top_wall2 = find_wall_indices(x, y, left2, right2, down2, up2)

# room 3
left3 = find_index(x, 10.) + 1
right3 = find_index(x, xmax)
down3 = find_index(y, 0.) + 1
up3 = find_index(y, ymax)

left_wall3, right_wall3, bottom_wall3, top_wall3 = find_wall_indices(x, y, left3, right3, down3, up3)

# grzejnik w room 1
h_left1 = 1
h_right1 = 2
h_down1 = 2
h_up1 = 4
area1 = (2-1) * (4-2)/hx**2

# grzejnik w room 2
h_left2 = 1
h_right2 = 2
h_down2 = 8
h_up2 = 10
area2 = (2-1) * (10-8)/hx**2

# grzejnik w room 2
h_left3 = 16
h_right3 = 17
h_down3 = 4
h_up3 = 6
area3 = (17-16) * (6-4)/hx**2

hpower = 2000
max_temp_heater = 20 + 273

grzejnik1 = Heater(h_left1, h_right1, h_down1, h_up1, hpower, max_temp_heater, area1)
grzejnik2 = Heater(h_left2, h_right2, h_down2, h_up2, hpower, max_temp_heater, area2)
grzejnik3 = Heater(h_left3, h_right3, h_down3, h_up3, hpower, max_temp_heater, area3)

pokoj1 = Room(left1, right1, down1, up1, 1, grzejnik1)
pokoj2 = Room(left2, right2, down2, up2, 1, grzejnik2)
pokoj3 = Room(left3, right3, down3, up3, 1, grzejnik3)

okno1 = Window(3, 6, 0, 0.5, 5+273)
okno2 = Window(0, 0.5, 8, 10, 5+273)
okno3 = Window(13, 16, 0, 0.5, 5+273)
okna = [okno1, okno2, okno3]

drzwi1 = Door(10, 10, 2.5, 3.5)
drzwi2 = Door(10, 10, 8.5, 9.5)
drzwi = [drzwi1, drzwi2]

neumann = [left_wall1, right_wall1, bottom_wall1, top_wall1, left_wall2, right_wall2, bottom_wall2, top_wall2,
           left_wall3, right_wall3, bottom_wall3, top_wall3]
grzejniki = [grzejnik1, grzejnik2, grzejnik3]
pokoje = [pokoj1, pokoj2, pokoj3]
dom = House(np.zeros((len(y), len(x))), pokoje)

alg = Solver(dom, ht, hx, x, y, X, Y, t, alpha, ro, c, 15+273, neumann, grzejniki, okna, drzwi)
alg.main_loop()
alg.plot_result()
alg.total_energy_use()
