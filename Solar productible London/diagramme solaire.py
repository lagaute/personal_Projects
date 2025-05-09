import numpy as np
import matplotlib.pyplot as plt

# Localisation de l'installation : LATITUDE PHI
phi = 51
phi = np.pi * phi / 180  # On passe tous les angles en radian pour les calculs avec les fonctions trigonométriques.

# Jour pour lequel on fait les calculs au format [ANNEE MOIS JOUR]
TIME = [2019,3, 20]

# On calcule le nombre de jours depuis le premier janvier à l'aide de la fonction datenum
jour = (np.datetime64(f'{TIME[0]}-{TIME[1]:02d}-{TIME[2]:02d}') - np.datetime64(f'{TIME[0]}-01-01')).astype(int) + 1


# On calcule la déclinaison en fonction du jour seulement
dec = 23.45 * np.sin(2 * np.pi * (jour + 284) / 365)  # en degré, equation très approximative
dec = np.pi * dec / 180  # en radian

def position_soleil(omega, delta, phi):
    h = np.arcsin(np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega))  # radian pour utilisation dans l'azimut
    a = (np.arcsin(np.cos(delta) * np.sin(omega) / np.cos(h))) * 180 / np.pi
    h = h * 180 / np.pi
    return h, a

# CALCUL SUR UNE JOURNEE
# On initialise un compteur qui s'incrémente à chaque passage dans la boucle
k = 0
hauteur = np.zeros(144)  # 24 heures * 6 intervals par heure = 144

for heure in range(24):
    for minutes in range(0, 60, 10):
        k = k + 1
        # Angle horaire omega.
        omega = -180 + (heure + minutes / 60) / 24 * 360
        omega = np.pi * omega / 180

        hauteur[k - 1], _ = position_soleil(omega, dec, phi)

t = np.linspace(0, 23 + 50 / 60, k)

plt.plot(t, hauteur)
plt.plot([0, 23], [0, 0], 'r')
plt.title("Hauteur du soleil en fonction de l'heure solaire locale")
plt.xlabel('Heure')
plt.ylabel('Hauteur en degré')
plt.show()