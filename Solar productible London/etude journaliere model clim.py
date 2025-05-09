import numpy as np
import matplotlib.pyplot as plt

# Localisation de l'installation : LATITUDE PHI
phi = 51

# Paramètres de l'installation Initialisation des valeurs propres au panneau
gamma_orientation = 10  # orientation du panneau (0° pour le sud, compté positivement vers l'ouest)
inclinaison = 0  # Inclinaison du panneau par rapport à l'horizontale

# On fixe le jour pour lequel on fait les calculs [ANNEE MOIS JOUR]
TIME = [2019, 3, 20]

albedo = 0.2  # Albédo

# Paramètres climatiques: Modèle Perrin de Brichambault
# ciel très pur
A = 1200
A_p = 187
A_p_p = 990
B = 2.5
B_p_p = 1.25
# ciel moyennement trouble
# A = 1230
# A_p = 125
# A_p_p = 1080
# B = 4
# B_p_p = 1.22
# ciel très trouble
# A = 1200
# A_p = 187
# A_p_p = 990
# B = 2.5
# B_p_p = 1.25

phi = np.pi * phi / 180  # On passe tous les angles en radian pour les calculs avec les fonctions trigonométriques.
gamma_orientation = np.pi * gamma_orientation / 180
inclinaison = np.pi * inclinaison / 180

# On calcule le nombre de jours depuis le premier janvier à l'aide de la
# fonction datenum
jour = (np.datetime64(f'{TIME[0]}-{TIME[1]:02d}-{TIME[2]:02d}') - np.datetime64(f'{TIME[0]}-01-01')).astype(int) + 1

# On calcule la déclinaison en fonction du jour seulement
declinaison = 23.45 * np.sin(2 * np.pi * (jour + 284) / 365)  # en degré, equation très approximative
declinaison = np.pi * declinaison / 180  # en radians


def irradiances(omega, delta, phi, gam, i, albedo, A, A_p, A_p_p, B, B_p_p):
    h = np.arcsin(np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega))
    a = np.arcsin(np.cos(delta) * np.sin(omega) / np.cos(h))

    # Calcul du rayonnement sur un capteur d'orientation et d'inclinaison donnes
    # Irradiance directe normale I* lorsque la hauteur du soleil est positive modélisé (Modèle Perrin de Brichambault)
    if h > 0:
        I_etoile = A * np.exp(-1. / (B * np.sin(h + 2 * np.pi / 180)))
    else:
        I_etoile = 0

    # Calcul de l'irradiance directe S* sur la surface inclinée du capteur à partir de I* modélisé
    if (np.cos(h) * np.sin(i) * np.cos(a - gam) + np.sin(h) * np.cos(i)) <= 0:  # soleil "derrière" le capteur
        S_etoile = 0
    else:
        S_etoile = I_etoile * (np.cos(h) * np.sin(i) * np.cos(a - gam) + np.sin(h) * np.cos(i))

    # Calcul de l'irradiance diffuse D* sur la surface inclinée du capteur à partir de D_h* et de G_h* modélisés
    if h > 0:
        D_etoile_horizontal = A_p * (np.sin(h) ** 0.4)  # irradiance horizontale diffuse Modèle P. De B.
        G_etoile_horizontal = A_p_p * (np.sin(h) ** B_p_p)  # irradiance horizontale globale Modèle P. De B.
        D_etoile = D_etoile_horizontal * ((1 + np.cos(i)) / 2) + albedo * G_etoile_horizontal * (
                    (1 - np.cos(i)) / 2)
    else:
        D_etoile = 0
        G_etoile_horizontal = 0

    return S_etoile, D_etoile, G_etoile_horizontal


# CALCUL SUR UNE JOURNEE
# On initialise un compteur qui s'incrémente à chaque passage dans la boucle
k = 0
S_etoile = np.zeros(144)  # 24 heures * 6 intervals par heure = 144
D_etoile = np.zeros(144)
G_etoile_horizontal = np.zeros(144)

for heure in range(24):
    for minutes in range(0, 60, 10):
        k = k + 1

        # Angle horaire omega.
        omega = -180 + (heure + minutes / 60) / 24 * 360  # en degrés
        omega = np.pi * omega / 180  # en radians

        S_etoile[k - 1], D_etoile[k - 1], G_etoile_horizontal[k - 1] = irradiances(
            omega, declinaison, phi, gamma_orientation, inclinaison, albedo, A, A_p, A_p_p, B, B_p_p)

# REMARQUE G_etoile_horizontal ne sert pas toujours
G_etoile = S_etoile + D_etoile
t = np.linspace(0, 23 + 50 / 60, k)

plt.figure()
plt.plot(t, G_etoile)
plt.title('Irradiance globale sur le capteur')
plt.xlabel('Heure')
plt.ylabel('Rayonnement global en W/m²')
plt.figure()
plt.plot(t, S_etoile)
plt.title('Irradiance directe sur le capteur')
plt.xlabel('Heure')
plt.ylabel('W/m²')
plt.figure()
plt.plot(t, D_etoile)
plt.title('Irradiance diffuse sur le capteur')
plt.xlabel('Heure')
plt.ylabel('Rayonnement diffus en W/m²')

# valeur Rg à midi
omega = 0
S_etoile_midi, D_etoile_midi, G_etoile_horizontal_midi = irradiances(0, declinaison, phi, gamma_orientation, inclinaison,
                                                                      albedo, A, A_p, A_p_p, B, B_p_p)
R_g_midi = (S_etoile_midi + D_etoile_midi) / G_etoile_horizontal_midi

plt.show()