#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NB_ATOM_NAMES 119
#define NB_MOLECULES 90130


struct nom_at { char c1, c2; };
struct liaison { int A1, A2; int l_type; };
struct isthme
{
	struct liaison l;
	int id_composant;
}isthme;
typedef struct isthme isthmes;

struct couplet{
	int a1;
	int a2;
	int poids;
	int type; //1 s'il s'agit d'une liaison distance et 0 si les deux cycles partagent des liaison
}couplet;

struct graphe{
	int nb_sommets;
	int nb_arete;
	int nb_connexe;
	struct lesommet *som;
	struct couplet *aretes;
	int **matrice_cycles_type;
	int **matrice_cycles_poids;
}graphe;
struct unsommet
{
	int id;
	int poids;
}unsommet;
typedef struct unsommet SOMMET;
struct graphemoleculaire
{
	int nb_atomes;
	int nb_liaisons;
	int *liste_atomes;
	int *type_atomes;
	struct liaison *liste_liaisons;
	int **matrice_liaisons;
	int nb_connexe;
	int pere;//son numero de generation : 0 graphe de la molécule initial 1 : composantes connexe du grape initial 2: sousgraphes obentus en retirant les isthmes
}graphemoleculaire;

typedef struct graphemoleculaire graphemol;
struct molecule{
	char chebi_name[1024];
	int chebi_id;
	int nb_atomes;
	int nb_hydrogene;
	int nb_liaisons;
	int *liste_atomes;
	int **matrice_liaisons;
	struct liaison *liste_liaisons;
	struct graphe g;
	int g_def;
} molecule;

struct uncycle
{
	int nb_atomes;
	int sommet;
	struct liaison arete;
	int *chemin1;
	int *chemin2;
	int pere;
	int *sommets;
	int id_cycle;
}uncycle;

typedef struct uncycle cycles;

int atom_num (char *name);
void init_atom_num ();
int lire_num_atome(FILE *F);
int valeur_char (FILE *F);
void ligne_suivante(FILE *F);
int lire_entier_3 (FILE * F);
struct liaison lire_liaison(FILE *F);
int lire_chebi_id(FILE *F);
void lire_chebi_name(FILE *F, struct molecule *M);
void lire_fin_molecule(FILE *F);
void trouver_la_fin_de_M(FILE *F);
struct molecule lire_molecule(FILE *F);
void tailles_molecules(struct molecule *M);
void liberer_molecule(struct molecule m);
