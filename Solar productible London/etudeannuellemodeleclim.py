import numpy as np
import matplotlib.pyplot as plt

# Calcul de l'irradiation solaire en un lieu donné.

# Début ENTREES UTILISATEURS
# Localisation de l'installation:LATUTUDE PHI
phi = 45

# Paramètres de l'installation Initialisation des valeurs propres au panneau
gamma_orientation = 0  # orientation du panneau (0° pour le sud, compté positivement vers l'ouest)
inclinaison = 10  # Inclinaison du panneau par rapport à l'horizontale

albedo = 0.2  # Albédo

# Paramètres climatiques: Modèle Perrin de Brichambault
# ciel très pur
A, A_p, A_p_p, B, B_p_p = 1230, 125, 1080, 4, 1.22
# ciel moyennement trouble
# A, A_p, A_p_p, B, B_p_p = 1230, 125, 1080, 4, 1.22
# ciel très trouble
# A, A_p, A_p_p, B, B_p_p = 1200, 187, 990, 2.5, 1.25

# Fin ENTREES UTILISATEURS

phi = np.pi * phi / 180  # On passe tous les angles en radian pour les calculs avec les fonctions trigonométriques.
gamma_orientation = np.pi * gamma_orientation / 180
inclinaison = np.pi * inclinaison / 180


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
    if (np.cos(h) * np.sin(i) * np.cos(a - gam) + np.sin(h) * np.cos(i)) <= 0:  # soleil "derriere" le capteur
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

    return S_etoile, D_etoile


S_annuel = 0
D_annuel = 0
G_annuel = 0
S_mensuel = np.zeros(12)
D_mensuel = np.zeros(12)
G_mensuel = np.zeros(12)
nombre_de_jours = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
jourdumois = np.array([15, 46, 75, 106, 136, 167, 197, 228, 259, 289, 310, 340])


for mois in range(12):
    S_journalier = np.zeros(nombre_de_jours[mois])
    D_journalier = np.zeros(nombre_de_jours[mois])
    G_journalier = np.zeros(nombre_de_jours[mois])

    for jour_mois in range(nombre_de_jours[mois]):
        # On fixe le jour pour lequel on fait les calculs
        jour = jourdumois[mois]
        # On calcule la déclinaison en fonction du jour seulement
        declinaison = 23.45 * np.sin(2 * np.pi * (jour + 284) / 365)  # en degré
        # equation très approximative, d'autres plus précises existent dans la littérature.
        # prise en compte possible de l'instant dans la déclinaison
        declinaison = np.pi * declinaison / 180  # en radian

        # On initialise un compteur qui s'incrémente à chaque passage dans la boucle
        k = 0
        for heure in range(24):
            S_horaire = np.zeros(24)
            D_horaire = np.zeros(24)
            G_horaire = np.zeros(24)

            for minutes in range(0, 60, 10):
                k += 1
                # Angle horaire omega.
                omega = -180 + (heure + minutes / 60) / 24 * 360
                omega = np.pi * omega / 180

                S_etoile_k, D_etoile_k = irradiances(omega, declinaison, phi, gamma_orientation, inclinaison, albedo,
                                                     A, A_p, A_p_p, B, B_p_p)

                G_etoile_k = S_etoile_k + D_etoile_k

                S_horaire[heure] += S_etoile_k / 6 / 1000  # Division par 6 car calcul par pas de 10min et par 1000 pour passage  en KWh
                D_horaire[heure] += D_etoile_k / 6 / 1000
                G_horaire[heure] += G_etoile_k / 6 / 1000

            S_journalier[jour_mois] += S_horaire[heure]
            D_journalier[jour_mois] += D_horaire[heure]
            G_journalier[jour_mois] += G_horaire[heure]

    # intégration mensuelle
    S_mensuel[mois] += np.sum(S_journalier)
    D_mensuel[mois] += np.sum(D_journalier)
    G_mensuel[mois] += np.sum(G_journalier)

S_annuel = np.sum(S_mensuel)  # KWh/m²
D_annuel = np.sum(D_mensuel)  # KWh/m²
G_annuel = np.sum(G_mensuel)  # KWh/m²

fig, axs = plt.subplots(3, 1, figsize=(10, 15))


axs[0].bar(range(1, 13), S_mensuel, color='r')
axs[0].set_title('Irradiation directe sur le capteur (KWh/m²)')
axs[0].set_xlabel('Mois')
axs[0].set_ylabel('KWh/m²')

axs[1].bar(range(1, 13), D_mensuel, color='b')
axs[1].set_title('Irradiation diffuse sur le capteur (KWh/m²)')
axs[1].set_xlabel('Mois')
axs[1].set_ylabel('KWh/m²')

axs[2].bar(range(1, 13), G_mensuel, color='g')
axs[2].set_title('Irradiation globale sur le capteur (KWh/m²)')
axs[2].set_xlabel('Mois')
axs[2].set_ylabel('KWh/m²')
plt.tight_layout()
plt.show()