import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as am


class Heater:
    def __init__(self, left, right, down, up, power, max_temp, area):
        self.left = left
        self.right = right
        self.down = down
        self.up = up
        self.power = power
        self.max_temp = max_temp
        self.area = area    # pole powierzchni grzejnika


class Window:
    def __init__(self, left, right, down, up, temperature):
        self.left = left
        self.right = right
        self.down = down
        self.up = up
        self.temperature = temperature  # temperatura na oknie


class Door:
    def __init__(self, left, right, down, up):
        self.left = left
        self.right = right
        self.down = down
        self.up = up


class Room:
    """ Reprezentuje pojedynczy pokój."""
    def __init__(self, left, right, down, up, u_prev, heater):
        self.left = left
        self.right = right
        self.down = down
        self.up = up
        self.u_prev = u_prev
        self.width = len(range(self.left, self.right))
        self.height = len(range(self.down, self.up))
        D2x = (np.diag(-2 * np.ones(self.height)) + np.diag(np.ones(self.height - 1), 1) +
               np.diag(np.ones(self.height - 1), -1))
        D2y = (np.diag(-2 * np.ones(self.width)) + np.diag(np.ones(self.width - 1), 1) +
               np.diag(np.ones(self.width - 1), -1))
        self.laplacian = np.kron(np.eye(self.width), D2x) + np.kron(D2y, np.eye(self.height))
        self.heater = heater
        self.energy_step = 0

    def update_u_prev(self, new_data):
        """Aktualizuje dane w pokoju"""
        if len(new_data) != len(self.u_prev):
            raise ValueError("Length of new data is not compatible!")
        self.u_prev = new_data

    def create_u_next(self, a, h_t, h_x, f):
        """ Wyznacza kolejny krok rozwiązania numerycznego dla pokoju. """
        termostat = 1
        if np.mean(self.u_prev) > self.heater.max_temp:      # sprawdza czy śr temp w pokoju <= tej ust. na grzejniku
            termostat = 0

        u_n = (self.u_prev + a * (h_t/h_x**2)*np.matmul(self.laplacian, self.u_prev) +
               h_t*termostat*f)
        self.energy_step = h_t*termostat*f.sum()
        return u_n


class House:
    """ Reprezentuje modelowane mieszkanie. """
    def __init__(self, big_matrix, rooms):
        self.big_matrix = big_matrix         # macierz reprezentująca plan domu
        self.rooms = rooms      # lista z instancjami 'room' z których składa się dom

    def divide_big_matrix(self):
        """Dzieli macierz domu na pokoje"""
        for r in self.rooms:
            r.update_u_prev(self.big_matrix[r.down:r.up, r.left:r.right].flatten())
            # 'poprawka' to -1 przy r.left i r.right powyżej

    def combine_into_big_matrix(self):
        """Łączy macierze pokoi w macierz dom"""
        for r in self.rooms:
            self.big_matrix[r.down:r.up, r.left:r.right] = r.u_prev.reshape(r.height, r.width)

    def solve_step(self, a, h_t, h_x, f):
        """Jeden krok przepływu ciepło w mieszkaniu"""
        self.divide_big_matrix()
        for r in self.rooms:
            f_r = f[r.down:r.up, r.left:r.right].flatten()
            r.update_u_prev(r.create_u_next(a, h_t, h_x, f_r,))
        self.combine_into_big_matrix()


class Solver:
    def __init__(self, house, h_t, h_x, x, y, X, Y, t, a, ro, c, temp0, neumann, heaters, windows, doors):
        self.house = house
        self.temp0 = temp0
        self.a = a
        self.ro = ro
        self.c = c
        self.h_t = h_t
        self.h_x = h_x
        self.x = x
        self.y = y
        self.X, self.Y = X, Y
        self.t = t
        self.u = np.zeros((len(self.x)*len(self.y), len(self.t)))
        self.heaters = heaters
        self.windows = windows
        self.doors = doors
        self.neumann = neumann
        self.psi = [0]        # do podliczania całkowitej energii

    def apply_dirichlet_conditions(self):
        """ Nakłada warunki Dirichleta na okna. """
        for window in self.windows:
            left_index = find_index(self.x, window.left)
            right_index = find_index(self.x, window.right)
            down_index = find_index(self.y, window.down)
            up_index = find_index(self.y, window.up)

            self.house.big_matrix[down_index:up_index, left_index:right_index] = window.temperature

    def connecting_doors(self):
        """ Uśrednia temperaturę na drzwiach"""
        for door in self.doors:
            left_index = find_index(self.x, door.left)
            right_index = find_index(self.x, door.right)
            down_index = find_index(self.y, door.down)
            up_index = find_index(self.y, door.up)

            self.house.big_matrix[down_index:up_index + 1, left_index:right_index + 1] = (self.house.big_matrix[
                                                                                          down_index:up_index + 1,
                                                                                          left_index - 1:right_index] +
                                                                                          self.house.big_matrix[
                                                                                          down_index:up_index + 1,
                                                                                          left_index + 1:right_index + 1 + 1]) / 2

    def apply_neumann_conditions(self):
        i = 0
        while i < len(self.neumann):
            # lewy brzeg
            self.house.big_matrix.flatten()[self.neumann[i]] = self.house.big_matrix.flatten()[(self.neumann[i] + 1)]
            # prawy brzeg
            self.house.big_matrix.flatten()[self.neumann[i + 1]] = self.house.big_matrix.flatten()[self.neumann[i + 1] - 1]
            # dolny brzeg
            self.house.big_matrix.flatten()[self.neumann[i + 2]] = self.house.big_matrix.flatten()[
                self.neumann[i + 2] + len(self.y)]
            # górny brzeg
            self.house.big_matrix.flatten()[self.neumann[i + 3]] = self.house.big_matrix.flatten()[
                self.neumann[i + 3] - len(self.y)]
            i += 4

    def prepare_force(self):
        f = np.zeros((len(self.y), len(self.x)))
        for h in self.heaters:
            X_heater, Y_heater = np.meshgrid(self.x, self.y)
            distance_squared = (X_heater - (h.left + h.right) / 2) ** 2 + (Y_heater - (h.down + h.up) / 2) ** 2
            sigma = 1
            A = 60 * h.power / h.area  # 60* bo czas w minutach?
            f += (A / (self.ro * self.c)) * np.exp(-distance_squared / (2 * sigma ** 2))   # * self.ro * self.c
            left_index = find_index(self.x, h.left)
            right_index = find_index(self.x, h.right)
            down_index = find_index(self.y, h.down)
            up_index = find_index(self.y, h.up)
            #f[down_index:up_index+1, left_index:right_index+1] = A/(self.ro * self.c)
        return f

    def main_loop(self):  # avg_temp_dist
        # outside_temp = avg_temp_dist[0]
        self.u[:, 0] = np.repeat(self.temp0, len(self.x)*len(self.y)).flatten()       # war. pocz., stała temp wszędzie
        force = self.prepare_force()
        for k, time in enumerate(self.t):
            if k == 0:
                for r in self.house.rooms:
                    r.u_prev = self.u[:, 0].reshape(len(self.y), len(self.x))[r.down:r.up, r.left:r.right].flatten()
                continue
            self.house.solve_step(self.a, self.h_t, self.h_x, force)

            # nakładanie warunków brzegowych
            self.apply_dirichlet_conditions()
            self.apply_neumann_conditions()

            # uśrednianie drzwi
            self.connecting_doors()

            # podliczenie energii w kroku
            psi_elem = 0
            for r in self.house.rooms:
                psi_elem += r.energy_step
            self.psi.append(psi_elem)

            self.u[:, k] = self.house.big_matrix.flatten()

    def plot_result(self):
        plt.figure(figsize=(10, 8))
        plt.pcolormesh(self.X, self.Y, self.u[:, -1].reshape(len(self.y), len(self.x)))
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Mapa ciepła")
        plt.colorbar()
        plt.show()

    def result_animated(self):
        fig, ax = plt.subplots(figsize=(10, 8))
        pcm = ax.pcolormesh(self.X, self.Y, self.u[:, -1].reshape(len(self.y), len(self.x)), shading='auto')
        colorbar = fig.colorbar(pcm, ax=ax)
        ax.set_title("Ewolucja mapy ciepła")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        # Funkcja aktualizująca ramki animacji
        def update(frame):
            data = self.u[:, frame].reshape(len(self.y), len(self.x))
            pcm.set_array(data.ravel())

            # Dynamiczna aktualizacja zakresów
            pcm.set_clim(vmin=data.min(), vmax=data.max())
            colorbar.update_normal(pcm)

            return pcm,

        step = 1
        frames = range(0+1, len(self.t), step)
        ani = am.FuncAnimation(fig, update, frames=frames, blit=True, interval=50)
        plt.show()

    def total_energy_use(self):
        total_energy = np.cumsum(self.psi) * self.h_x**2
        plt.plot(self.t, total_energy)
        plt.xlabel("czas")
        plt.ylabel("energia [J]")
        plt.title("Całkowita energia grzania")
        plt.show()


def find_index(arr, value):
    idx = np.where(arr == value)[0]
    if len(idx) > 0:
        return int(idx[0])
    else:
        raise ValueError(f"Nie znaleziono wartości {value} w tabeli.")