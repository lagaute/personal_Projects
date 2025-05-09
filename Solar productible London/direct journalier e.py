import numpy as np
import matplotlib.pyplot as plt

# Localisation de l'installation : LATITUDE PHI
phi = 51

# Paramètres de l'installation Initialisation des valeurs propres au panneau
gamma_orientation = 10  # orientation du panneau (hémisphère nord: 0° pour le sud, compté positivement vers l'ouest/hémisphère sud : 0° vers le nord)
inclinaison = 0  # Inclinaison du panneau par rapport à l'horizontale

# On fixe le jour pour lequel on fait les calculs [ANNEE MOIS JOUR]
TIME = [2019, 12, 21]

phi = np.pi * phi / 180  # On passe tous les angles en radian pour les calculs avec les fonctions trigonométriques.
gamma_orientation = np.pi * gamma_orientation / 180
inclinaison = np.pi * inclinaison / 180

# On calcule le nombre de jours depuis le premier janvier à l'aide de la fonction timedelta
jour = (np.datetime64(f'{TIME[0]}-{TIME[1]:02d}-{TIME[2]:02d}') - np.datetime64(f'{TIME[0]}-01-01')).astype(int) + 1

# On calcule la déclinaison en fonction du jour seulement
declinaison = 23.45 * np.sin(2 * np.pi * (jour + 284) / 365)  # en degré, equation très approximative
declinaison = np.pi * declinaison / 180  # en radians


def irradiances(omega, delta, phi, gam, i):
    h = np.arcsin(np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega))
    a = np.arcsin(np.cos(delta) * np.sin(omega) / np.cos(h))

    # Calcul du rayonnement sur un capteur d'orientation et d'inclinaison donnés
    # Irradiance directe lorsque la hauteur du soleil est positive (soleil au-dessus de l'horizon!)
    if h > 0:
        I_etoile = 1000
    else:
        I_etoile = 0

    # Calcul du flux (S) sur la surface inclinée du capteur
    if (np.cos(h) * np.sin(i) * np.cos(a - gam) + np.sin(h) * np.cos(i)) <= 0:  # soleil "derrière" le capteur
        S_etoile = 0
    else:
        S_etoile = I_etoile * (np.cos(h) * np.sin(i) * np.cos(a - gam) + np.sin(h) * np.cos(i))

    return S_etoile


# CALCUL SUR UNE JOURNEE
# On initialise un compteur qui s'incrémente à chaque passage dans la boucle
k = 0
S_etoile_journalier = np.zeros(144)  # 24 heures * 6 intervals par heure = 144

for heure in range(24):
    for minutes in range(0, 60, 10):
        k = k + 1

        # Angle horaire omega.
        omega = -180 + (heure + minutes / 60) / 24 * 360  # en degrés
        omega = np.pi * omega / 180  # en radians

        S_etoile_journalier[k - 1] = irradiances(omega, declinaison, phi, gamma_orientation, inclinaison)

t = np.linspace(0, 23 + 50 / 60, k)

plt.plot(t, S_etoile_journalier)
plt.title("Irradiance du capteur en fonction de l'heure solaire locale")
plt.xlabel('Heure')
plt.ylabel('Irradiance en W/m^2')
plt.show()