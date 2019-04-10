#include "modelisation.h"
#include "string.h"
#define taille 10000
#define atomes 10000

int last_chrono;

int main(int argc, char *argv[])
{
	if( argc != 2 )
	{
		fprintf(stdout,"Missing arguments (CHEBI ID) et taille max des cycles \n");
		exit(20);
	}
	//on recupère le chebi de la molécule passée en parametre

	taille_cycle = atoi(argv[1]);
	

	init_atom_num ();
	struct molecule *M = lecture_fichier_chebi();
	


	

	char src[64];
	sprintf(src,"stats/tailles_cycles.data");
	
	FILE *F = fopen(src,"w");

	if( F == NULL)
	{
		fprintf(stdout,"Cannot create the file %s\n",src);
		exit(87);
	}
	
	sprintf(src,"stats/moyenne_cycles.data");
	FILE *G = fopen(src,"w");

	if( G == NULL)
	{
		fprintf(stdout,"Cannot create the file %s\n",src);
		exit(88);
	}
	
	sprintf(src,"stats/distribution_cycles.data");
	FILE *H = fopen(src,"w");

	if( H == NULL)
	{
		fprintf(stdout,"Cannot create the file %s\n",src);
		exit(89);
	}
	sprintf(src,"stats/nombre_atomes.data");
	FILE *L = fopen(src,"w");
	if( L == NULL)
	{
		fprintf(stdout,"Cannot create the file %s\n",src);
		exit(90);
	}
	
	sprintf(src,"stats/distribution_atomes.data");
	FILE *N = fopen(src,"w");

	if( N == NULL)
	{
		fprintf(stdout,"Cannot create the file %s\n",src);
		exit(91);
	}
	
	sprintf(src,"stats/moyenne_atomes.data");
	FILE *O = fopen(src,"w");

	if( O == NULL)
	{
		fprintf(stdout,"Cannot create the file %s\n",src);
		exit(91);
	}
	
	
	
	int i,somme = 0;
	int tab[taille];
	for( i = 0; i < taille;i++)
	{
		tab[i] = 0 ;
	}
	
	int tableau[atomes], moyenne = 0;
	for( i = 0; i < atomes;i++)
	{
		tableau[i] = 0 ;
	}
	
	
	for( i = 0; i < NB_MOLECULES; i++)
	{
		
		if (i % 50 == 0) 
		{ 
			fprintf(stdout,"\r%5d / %d",i,NB_MOLECULES);
			fflush(stdout); 
		}
		nb_arete_base = 0;
		taille_base = 0;
		GRAPHE_CYCLE c = construction_graphe_cycles(M[i]);
		fprintf(F,"%d %d\n",M[i].chebi_id, c.nb_sommets);
		fprintf(L,"%d %d\n",M[i].chebi_id, M[i].nb_atomes);
		somme +=c.nb_sommets;
		tab[c.nb_sommets]++;
		fprintf(L,"%d %d\n",M[i].chebi_id, M[i].nb_atomes);
		moyenne += M[i].nb_atomes;
		tableau[M[i].nb_atomes]++;
		liberer_graphe_cycles(c);
	}
	
	int mols = 0;
	for( i = 0; i < taille; i++)
	{
		mols += tab[i];
		fprintf(H,"%d %d\n",i,tab[i]);
		if(mols == NB_MOLECULES)
			break;
		
	}
	int ats = 0;
	for( i = 0; i < atomes;i++)
	{
		ats += tableau[i];
		fprintf(N,"%d %d\n",i,tableau[i]);
		if(ats == NB_MOLECULES)
			break;
	}
	
	fprintf(O,"%f", (float) moyenne / NB_MOLECULES);
	
	fprintf(G,"%f\n",(float)somme/NB_MOLECULES);
	
	
	fclose(F);
	fclose(G);
	fclose(H);
	fclose(L);
	fclose(N);
	fclose(O);
	
	
	int nb_mol;
	printf("\n3. Libération de la mémoire : %.3lf s\n",chrono());

	for(nb_mol=0 ; nb_mol < NB_MOLECULES ; nb_mol++) 
		liberer_molecule(M[nb_mol]);
	free(M);

	system("gnuplot distribution.gplt");

	exit(0);
}
