#include "lecture_molecule_sdf.h"

struct couple{
	int a1;
	int a2;
}couple;

struct liste_voisins{
	int id_atome;
	int nb_voisins;
	int *id_voisins;
		
}liste_voisins;

struct distance{ 
	
	int *sommets;
	int *predecesseur;
}distance;

struct cycle{
	int poids;
	int *c;
}cycle;

struct vecteur{
	int taille;
	struct couple *sommets;
}vecteur;

struct lesommet{
	int id;
	int taille;
}lesommet;
struct arete_base
{
	int id1;
	int id2;
	int type;
	int poids;
}arete_base;
typedef struct arete_base ARETE;

struct graphe_cycle
{
	int nb_sommets;
	int nb_aretes;
	ARETE *liste_aretes;
	SOMMET *liste_sommets;
}graphe_cycle;
typedef struct graphe_cycle GRAPHE_CYCLE;

char *atom_name[NB_ATOM_NAMES];
int chebi_id;
int taille_cycle;
int nb_arete_base;
ARETE *base_aretes;
int distance_LC;
int *arete_cycle;
int **arete_liste;
int taille_base;
cycles *labase;
//liste des fonctions de : fonctions_modelisation.c
double chrono();
graphemol suppression_aretes_mol(graphemol m, int i,int *degre);
int verif_mol(graphemol g, int* degre, int *deja_elimine);
graphemol modification_structure_mol(graphemol g, int* degre,int* deja_elimine);
int *calcul_degre_mol( graphemol m);
void affichage_degre(struct molecule m, int *degre);
void affiche_graphemol(graphemol g);
graphemol conversion_mol_graphe(struct molecule m);
void liberer_graphemol(graphemol g);
void trouver_squelette_cycles();
void elimination_feuilles(struct molecule m);
struct molecule construction_matrice_mol(struct molecule m);
void affiche_matrice(struct molecule m);
void affiche_mol(struct molecule M);
void trouver_squelette_cycles(graphemol m);
void liberer_memoire_voisins(struct liste_voisins *v,graphemol m);
int *plus_court_chemin(int sommet1,int sommet2,graphemol m);
void affiche_pcc_chemin(int sommet1,int sommet2,int *chemin ,graphemol m);
int chemin_independant(int *chemin1, int *chemin2,struct liaison l);
int verification_ajout_cycle(cycles *l, int nb_cycles , cycles c,graphemol m);
cycles creer_un_cycle(graphemol m , int sommet, struct liaison l, int *chemin1,int *chemin2);
void afficher_un_cycle(cycles c, graphemol m );
void liberer_un_cycle(cycles c);
cycles *ajouter_un_cycle(cycles *liste, int nb_cycles , cycles c,graphemol m);
void obtenir_la_base(cycles *liste,int nb_cycles, graphemol m);
int position_de_arete(int sommet1, int sommet2,graphemol m);
int *concatener_deux_chemins(cycles c);
void affiche_matrice_cycles ( int **matrice , int nb_l,int nb_c);
int verification_egalite_tableaux(int *tab1,int *tab2 , int taille);
int fonction_xor(int a , int b);
int *produit_xor_matrice(int *t, int *tab1, int *tab2,int nb);
void ajouter_cycle_base(cycles c,graphemol m);
void fichier_base();
void afficher_un_cycle_fichier(cycles c,FILE *f);
int sommet_dans_basniveau(graphemol t,int p);
int distance_inter_magma(cycles a , cycles b , struct molecule m);
void nouvelle_arete(ARETE a);
int sommets_commun ( cycles a, cycles b);
int calcul_LC(cycles a, cycles b,struct molecule m);
int nombre_de_cycles(int sommet1,int *chemin,int sommet2,graphemol g);
void arete_dans_cycle(struct molecule m);
int position_graphemol_arete(graphemol g, int sommet1,int sommet2);
void arete_dans_cycle_liste(struct molecule m);
void fichier_dot();
int verification_LC(cycles a, cycles b,struct molecule m,int dist);
GRAPHE_CYCLE construction_graphe_cycles(struct molecule m);
void liberer_graphe_cycles( GRAPHE_CYCLE c);
ARETE copier_arete(ARETE a , ARETE b);


float similarite(GRAPHE_CYCLE a,GRAPHE_CYCLE b);
GRAPHE_CYCLE ajouter_un_sommet(int a , int b, int taille ,GRAPHE_CYCLE c);
GRAPHE_CYCLE construction_graphe_produit(GRAPHE_CYCLE a, GRAPHE_CYCLE b);
int min( int a , int b);
int valeur_absolue(int a);
int poids_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b);
int type_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b);
GRAPHE_CYCLE ajouter_une_arete_graphe(int a , int b, int arete ,GRAPHE_CYCLE c);
int ** construction_matrice_produit(GRAPHE_CYCLE p);
void liberation_matrice(int **matrice, int taille);
void affiche_matrice_produit(int **matrice, int taille);
void la_clique_max( int **matrice, int sommets,int **depart,GRAPHE_CYCLE produit);
void calcul_clique(int **matrice,int sommets,int *dans_clique,int taille_clique,int *candidat,int taille_candidat, double date,GRAPHE_CYCLE produit, int **depart);
int somme_sommets_aretes( int *clique, GRAPHE_CYCLE produit, int **depart, int sommets);
//liste des fonctions dans connexite.c

struct liste_voisins*  construction_voisinage_graphemol(graphemol M);
void affichage_liste_voisinage_molecule(struct liste_voisins* voisins,struct molecule M);
void affichage_liste_voisinage_graphemol(struct liste_voisins* voisins,graphemol M);
int position_graphemol(graphemol g,int element);
int existe_chaine_graphemol(int i, int j , graphemol m,struct liste_voisins * v);
int est_connexe_graphemol(graphemol m, int *deja_elimine);
int * sommets_atteignable_graphemol(int i,graphemol m,struct liste_voisins * v,int *deja_elimine);
void liberation_liste_voisins_graphemol(struct liste_voisins *v,graphemol m);
void liberation_graphe(struct graphe g);
void probleme_memoire();
graphemol * ensemble_connexe_graphemol(graphemol m , int *deja_elimine);
graphemol enlever_une_arete(graphemol g, struct liaison l);
graphemol ajouter_une_arete(graphemol g, struct liaison l);
int retirer_liaison_connexe(graphemol g, struct liaison l, int *deja_elimine);
int nombre_isthmes_graphemol(graphemol g, int *deja_elimine);
isthmes *retrouver_tous_isthmes(graphemol *liste_connexe,isthmes *liste, int *deja_elimine);
void affiche_liste_isthmes ( isthmes *l, int nb);
graphemol enlever_tous_isthmes(graphemol g,int *deja_elimine,isthmes *liste_isthmes,int nb_isthmes);
struct molecule * lecture_fichier_chebi();
struct molecule lire_molecule_sdf(FILE *F);
int position_M( int g1_chebi,struct molecule *M);
























