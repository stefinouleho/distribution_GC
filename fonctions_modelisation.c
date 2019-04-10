#include "modelisation.h"
#include <sys/time.h>
#define AUCUNE_LIAISON (-1024)


int taille_clique_max =0;
int *dans_clique_max= NULL;
int liaison_max = 0;
double last_chrono;
double date_max = 100000;

double chrono() 
{
	struct timeval tv;
	static double date_deb = 0.0;
	double date_courante;
	gettimeofday(&tv,NULL);
	date_courante = tv.tv_sec + ((double)tv.tv_usec)*1e-6;
	if (date_deb == 0.0) date_deb = date_courante;
	return date_courante-date_deb;
}



graphemol suppression_aretes_mol(graphemol m, int i,int *degre)
{
	graphemol M;
	if(degre[i - 1] == 1)
		M.nb_liaisons = m.nb_liaisons - 1;
	else
		M.nb_liaisons = m.nb_liaisons;
	M.nb_atomes = m.nb_atomes;
	M.liste_atomes   = malloc(M.nb_atomes   * sizeof(int));
	M.type_atomes   = malloc(M.nb_atomes   * sizeof(int));
	M.liste_liaisons = malloc(M.nb_liaisons * sizeof(struct liaison));
	M.matrice_liaisons = malloc(M.nb_atomes   * sizeof(int *));
	M.nb_connexe = m.nb_connexe;
	M.pere = m.pere;
	int t;
	for ( t = 0; t < m.nb_atomes;t++)
		M.matrice_liaisons[t] = malloc(M.nb_atomes *sizeof(int));
	int a,k = 0; 
	//mise a jour des atomes
	for (a = 0; a < m.nb_atomes ; a++) 
	{
		
			M.type_atomes[k] = m.type_atomes[a];
			M.liste_atomes[k] = m.liste_atomes[a];
			k++;
	}	
	for ( t = 0; t < m.nb_atomes;t++)
	{
		for ( a = 0; a < M.nb_atomes; a++)
			M.matrice_liaisons[t][a] = 0;
	}
	k = 0;
	for (a = 0; a < m.nb_liaisons; a++) 
	{
		if(m.liste_liaisons[a].A1 != i && m.liste_liaisons[a].A2 != i )
		{
			M.liste_liaisons[k].A1 = m.liste_liaisons[a].A1;
			M.liste_liaisons[k].A2 = m.liste_liaisons[a].A2 ;
			M.liste_liaisons[k].l_type = m.liste_liaisons[a].l_type;
			M.matrice_liaisons[M.liste_liaisons[k].A1 - 1][M.liste_liaisons[k].A2 -1] = M.liste_liaisons[k].l_type;
			k++;
		}
	}	
	liberer_graphemol(m);
	return M;
}

int verif_mol(graphemol g, int* degre, int *deja_elimine)
{
	int i ,verif = -1; 

	for ( i = 0; i < g.nb_atomes; i++)
	{
		if (degre[i] == 1 || degre[i] == 0 )
		{
			if(deja_elimine[i] == 0)
			{
				verif = i;
				deja_elimine[i] = 1;
				break;
			}
			
		}
			
	}
	return verif;
}

graphemol modification_structure_mol(graphemol g, int* degre,int* deja_elimine)
{
	int i = verif_mol(g,degre,deja_elimine) ;
	//affichage_degre(m,degre); 
	while(i != -1)
	{
		//fprintf(stdout, "le sommet %d\n",i +1  );
		g = suppression_aretes_mol( g, i + 1,degre);

		if ( degre != NULL)
			free(degre);
		degre = calcul_degre_mol(g);
		i = verif_mol(g,degre,deja_elimine) ; 
	}
	free(degre);
	return g;
}

int *calcul_degre_mol( graphemol m)

{
	int *degre = malloc( m.nb_atomes * sizeof( int)); // le degré de tous les sommets du graphe

	int i; 
	for(i = 0; i < m.nb_atomes; i++)
		degre[i] = 0;
	//remplissage
	int pos1,pos2;
	for ( i = 0; i < m.nb_liaisons; i++)
	{
		pos1 = position_graphemol(m,m.liste_liaisons[i].A1);
		pos2 = position_graphemol(m,m.liste_liaisons[i].A2);
		degre[pos1]++;
		degre[pos2]++;
	}

	return degre;
}
void affichage_degre(struct molecule m, int *degre)
{
	fprintf(stdout, "Debut affichage liste des degrés\n");
	int i;
	for( i = 0; i < m.nb_atomes; i++)
		printf("%d %d\n",i+1, degre[i]);
	fprintf(stdout, "fin affichage liste des degrés\n" );
}
void affiche_graphemol(graphemol g)
{
	
	printf("---Affichage du graphe de la molécule ------ \n nb atomes %d \n nb liaisons %d \n", g.nb_atomes,g.nb_liaisons);
	printf("Generation : %d\n", g.pere);
	int i,j;
	for (i = 0; i < g.nb_atomes; i++)
	{
		printf("%2d ", g.liste_atomes[i]);
	}
	printf("\n");
	for (i = 0; i < g.nb_atomes; i++)
	{
		printf("%2d ", g.type_atomes[i]);
	}
	printf("\nMatrice de liaison \n");
	for (i = 0; i < g.nb_atomes; i++)
	{
		for (j = 0; j < g.nb_atomes; j++)
		{
			printf("%2d ", g.matrice_liaisons[i][j]);
		}
		printf("\n");
	}
	printf("Fin de l'affichage \n");
}

graphemol conversion_mol_graphe(struct molecule m)
{
	graphemol g;
	g.nb_atomes = m.nb_atomes;
	g.nb_liaisons = m.nb_liaisons;
	g.liste_atomes = malloc(g.nb_atomes *sizeof(int));
	g.type_atomes = malloc(g.nb_atomes *sizeof(int));
	g.liste_liaisons = malloc(g.nb_liaisons *sizeof(struct liaison));
	g.nb_connexe = 0;
	g.pere = 0;
	if(g.liste_atomes == NULL)
		probleme_memoire();
	if(g.type_atomes == NULL)
		probleme_memoire();
	
	int i,j;
	for (i = 0; i < g.nb_atomes; i++)
	{
		g.liste_atomes[i] = i + 1;
		g.type_atomes[i] = m.liste_atomes[i];
	}
	for (i = 0; i < g.nb_liaisons; i++)
	{
		g.liste_liaisons[i].A1 = m.liste_liaisons[i].A1;
		g.liste_liaisons[i].A2 = m.liste_liaisons[i].A2;
		g.liste_liaisons[i].l_type = m.liste_liaisons[i].l_type;
	}
	g.matrice_liaisons = malloc(g.nb_atomes *sizeof(int *));
	if( g.matrice_liaisons == NULL)
		probleme_memoire();
	for (i = 0; i < g.nb_atomes; i++)
	{
		g.matrice_liaisons[i]= malloc(g.nb_atomes *sizeof(int));
		if(g.matrice_liaisons[i] == NULL)
			probleme_memoire();
	}

	for (i = 0; i < g.nb_atomes; i++)
	{
		for (j = 0; j < g.nb_atomes; j++)
		{
			g.matrice_liaisons[i][j] = m.matrice_liaisons[i][j];
		}
	}
	//affiche_graphemol(g);
	return g;
}

void liberer_graphemol(graphemol g)
{
	
	if ( g.liste_atomes != NULL)
		free(g.liste_atomes);
	if ( g.type_atomes != NULL)
		free(g.type_atomes);
	if( g.liste_liaisons != NULL)
		free(g.liste_liaisons);
	
	int i;
	for( i = 0; i < g.nb_atomes;i++)
	{
		if(g.matrice_liaisons[i] != NULL)
			free(g.matrice_liaisons[i]);
	}
	if(g.matrice_liaisons != NULL)
		free(g.matrice_liaisons);
	

}
//liberer la memoire de liste voisins 

void liberer_memoire_voisins(struct liste_voisins *v,graphemol m)
{
	int i;
	for (i = 0; i < m.nb_atomes; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);

}
//calcule un plus court chemin entre deux sommets 

int *plus_court_chemin(int sommet1,int sommet2,graphemol m)
{
	int *chemin =  NULL;
	struct liste_voisins *v = construction_voisinage_graphemol(m);
	//affichage_liste_voisinage_graphemol(voisins,m);
	int distance[m.nb_atomes];
	int predecesseur[m.nb_atomes];

	int i,j;
	for (i = 0; i < m.nb_atomes; i++)
	{
		distance[i] = -1;
		predecesseur[i] = -1;
	}
		
	distance[position_graphemol(m,sommet1)] = 0;
	predecesseur[position_graphemol(m,sommet1)] = position_graphemol(m,sommet1);

	int pos = position_graphemol(m,sommet2);
	int iteration = 0;
	while(distance[pos] == -1)
	{
		for( i = 0; i < m.nb_atomes; i++)
		{
			if(distance[i] == iteration)
			{
				for(j = 0; j < v[i].nb_voisins;j++)
				{
					if(distance[position_graphemol(m,v[i].id_voisins[j])] == -1)
					{
						distance[position_graphemol(m,v[i].id_voisins[j])] = iteration + 1;
						predecesseur[position_graphemol(m,v[i].id_voisins[j])] = m.liste_atomes[i];
					}
				}
			}
		}
		iteration++;
	}
	chemin = malloc((distance[pos] +1) * sizeof(int));
	int s  = 0;
	chemin[s] = distance[pos];
	int pred = predecesseur[pos];
	s++;
	for(i = 0; i < distance[pos] - 1;i++)
	{
		chemin[s] = pred;
		pred  = predecesseur[position_graphemol(m,pred)];
		s++;
	}
	liberer_memoire_voisins(v,m);
	return chemin;
}

//affihce un plus court chemin
void affiche_pcc_chemin(int sommet1,int sommet2,int *chemin ,graphemol m)
{
	int i;
	printf("affiche un plus court chemin entre %d et %d  \n",sommet1,sommet2);
	printf("affiche un plus court chemin entre %d et %d de taille : %d \n",sommet1,sommet2 ,chemin[0] - 1);
	printf("%d ",sommet1);
	for( i  = chemin[0] - 1; i > 0; i--)
	{
		printf("%d ",chemin[i]);
	}
	printf("%d ",sommet2);
	printf("\n");

}

// retourne 1 si les chemins sont independants et 0 sinon 
int chemin_independant(int *chemin1, int *chemin2,struct liaison l)
{
	
	int resultat = 1,trouve = 0;
	int i,j;
	for( i = 0; i < chemin1[0] -1;i++)
	{
		if(chemin1[i+1] == l.A1 ||chemin1[i+1] == l.A2)
		{
			resultat = 0;
			return resultat;
		}
	}

	for( i = 0; i < chemin2[0] -1;i++)
	{
		if(chemin2[i+1] == l.A1 || chemin2[i+1] == l.A2)
		{
			resultat = 0;
			return resultat;
		}
	}
	if(chemin1[0] > chemin2[0])
	{

		for( i = 0 ; i < chemin2[0] - 1; i++)
		{
			for( j = 0; j < chemin1[0] - 1; j++)
			{
				if(chemin1[j + 1] == chemin2[i + 1])
				{
					resultat = 0;
					trouve = 1;
					break;
				}
			}
			if( trouve  == 1)
				break;
		}
	}
	else
	{
		for( i = 0 ; i < chemin1[0] - 1; i++)
		{
			for( j = 0; j < chemin2[0] - 1; j++)
			{
				if(chemin1[i + 1] == chemin2[j + 1])
				{
					resultat = 0;
					trouve = 1;
					break;
				}
			}
			if( trouve  == 1)
				break;
		}

	}
	return resultat;
}
//verifie si le cycle n'a pas encore été trouvé. retourne 1 si le cycle est deja das la base et 0 sinon
int verification_ajout_cycle(cycles *l, int nb_cycles , cycles c,graphemol m )
{
	int trouve = 0,i,j;
	
	int compteur ;	
	for( i =  0; i < nb_cycles; i++)
	{
		
		if(l[i].nb_atomes == c.nb_atomes)
		{
			compteur = 0;
			for( j =  0; j < c.nb_atomes; j++)
			{
				if(c.sommets[j] == l[i].sommets[j])
				{
					compteur++;
				}
			}
			if( compteur == c.nb_atomes)
			{
				trouve  = 1;
				break;
			}
		}
	}
	return trouve;
}
//cree le cycle qui passe le sommet et l'arete dans le graphe m
cycles creer_un_cycle(graphemol m , int sommet, struct liaison l, int *chemin1,int *chemin2)
{
	cycles c;
	int tableau[m.nb_atomes],i;

	for (i = 0; i < m.nb_atomes; i++)
		tableau[i] = 0;

	tableau[position_graphemol(m,sommet)] = 1;
	tableau[position_graphemol(m,l.A1)] = 1;
	tableau[position_graphemol(m,l.A2)] = 1;
	for(i  = 0; i < chemin1[0] - 1; i++)
		tableau[position_graphemol(m, chemin1[i+1])] = 1;
	for(i  = 0; i < chemin2[0] - 1; i++)
		tableau[position_graphemol(m, chemin2[i+1])] = 1;

	c.nb_atomes = 0;

	for (i = 0; i < m.nb_atomes; i++)
		c.nb_atomes+= tableau[i] ;

	c.sommet = sommet;
	c.arete  = l;
	c.chemin1 = malloc(chemin1[0] * sizeof(int));
	for( i = 0; i < chemin1[0]; i++)
		c.chemin1[i] = chemin1[i];
	
	c.chemin2 = malloc(chemin2[0] * sizeof(int));
	for( i = 0; i < chemin2[0]; i++)
		c.chemin2[i] = chemin2[i];
	
	c.sommets = malloc(c.nb_atomes * sizeof(int));
	int k = 0;
	for ( i = 0; i < m.nb_atomes; i++)
	{
		if(tableau[i] == 1)
		{
			c.sommets[k] = m.liste_atomes[i];
			k++;
		}
	}
	c.pere = 0;
	c.id_cycle = 0;
	return c;
}

//liberer un cycle

void liberer_un_cycle(cycles c)
{

	
	if( c.sommets != NULL)
		free(c.sommets);
	if(c.chemin1 != NULL)
		free(c.chemin1);
	if(c.chemin2 != NULL)
		free(c.chemin2);
}
//afficher un cycle 

void afficher_un_cycle(cycles c, graphemol m )
{
	printf("afficher le cycle\n");
	int i;
	printf("id_cycle : %d et pere = %d nb atomes %d , sommet = %d et liaison = [%d - %d]\n",c.id_cycle,c.pere,c.nb_atomes,c.sommet,c.arete.A1,c.arete.A2 );
	printf("Le cycle : %s%d -",atom_name[m.type_atomes[position_graphemol(m,c.sommet)]],c.sommet);
	for(i  = c.chemin1[0] - 1; i > 0; i--)
		printf(" %s%d -",atom_name[m.type_atomes[position_graphemol(m,c.chemin1[i])]],c.chemin1[i]);
	printf(" %s%d  - %s%d -", atom_name[m.type_atomes[position_graphemol(m,c.arete.A1)]],c.arete.A1,atom_name[m.type_atomes[position_graphemol(m,c.arete.A2)]],c.arete.A2);
	for(i  = 0; i < c.chemin2[0] - 1; i++)
		printf(" %s%d -",atom_name[m.type_atomes[position_graphemol(m,c.chemin2[i+1])]],c.chemin2[i+1]);
	printf(" %s%d\n",atom_name[m.type_atomes[position_graphemol(m,c.sommet)]],c.sommet);
}
//construction du graphe de cycles
void trouver_squelette_cycles(graphemol m)
{
	cycles *liste = NULL;
	cycles c;
	int nb_cycles  = 0; 
	int i,j;
	struct liaison l;
	int *chemin1,*chemin2;
	//affiche_graphemol(m);
	for (i = 0; i < m.nb_atomes; i++)
	{
		for(j = 0; j < m.nb_liaisons;j++)
		{
			l.A1 = m.liste_liaisons[0].A1;
			l.A2 = m.liste_liaisons[0].A2;
			l.l_type = m.liste_liaisons[0].l_type;
			m = enlever_une_arete(m,m.liste_liaisons[0]);
			//affiche_graphemol(m);	

			if(( l.A1 != m.liste_atomes[i])&&(l.A2 !=m.liste_atomes[i]))
			{
				//affiche_graphemol(m);
				chemin1 = plus_court_chemin(m.liste_atomes[i],l.A1,m);
				chemin2 = plus_court_chemin(m.liste_atomes[i],l.A2,m);
				if(chemin_independant(chemin1,chemin2,l))
				{
					c = creer_un_cycle(m,m.liste_atomes[i],l,chemin1,chemin2);
					c.pere = m.nb_connexe;
					if(!verification_ajout_cycle(liste,nb_cycles, c,m))
					{
						liste = ajouter_un_cycle(liste,nb_cycles,c,m);
						nb_cycles++;
					}
					liberer_un_cycle(c);
					
				}
					
				free(chemin1);
				free(chemin2);
				
			}
			m = ajouter_une_arete(m,l);
		}
	}

	obtenir_la_base(liste,nb_cycles,m);
	if( liste != NULL)
	{
		for ( i = 0; i < nb_cycles; i++)
		{
			liberer_un_cycle(liste[i]);
		}
			
		free(liste);
	}
		


}
void fichier_dot()
{

	char src[64];
	sprintf(src,"../CHEBI/%d/%d_squelette_LC_%d_%d.dot",chebi_id,chebi_id,taille_cycle,distance_LC);
	if(taille_base > 0)
	{
		FILE *D = fopen(src, "w");

		if ( D == NULL)
		{
			fprintf(stderr, "Cannot create the file %s \n", src);
			exit(45);
		}

		fprintf(D, "graph G {\n");
		fprintf(D, "node [fixedsize = true,width = 0.2,height = 0.2,fontsize = 5, shape=circle];\n edge [fontsize = 8];\n");
			
		int i;
		for ( i = 0 ; i < taille_base ; i++)
		{

			fprintf(D, "%d [ label=\"%d\" ];\n",labase[i].id_cycle,labase[i].nb_atomes);
		}

		for ( i = 0 ; i < nb_arete_base ; i++)
		{
			fprintf(D,"%d -- %d",base_aretes[i].id1,base_aretes[i].id2);
				
			if ( base_aretes[i].type == 1)
			{
				fprintf(D,"[color = blue,label = \"%d\"];\n",base_aretes[i].poids);
			}
			else if(base_aretes[i].type == 2)
			{
				fprintf(D,"[color = red,label = \"%d\"];\n",base_aretes[i].poids);
			}
			else
			{
				fprintf(D,"[color = green,label = \"%d\"];\n",base_aretes[i].poids);
			}
		}
				
			
		fprintf(D,"}");
		fclose(D);
	}	
	
}
void obtenir_la_base( cycles *liste,int nb_cycles, graphemol m)
{
	
	// classer les cyxcles par taille de poids ( base miniale selon horton)
	int i,j,l;

	int position[nb_cycles];
	int newposition[nb_cycles];
	for ( i = 0; i < nb_cycles;i++)
	{
		position[i] = liste[i].nb_atomes;
		newposition[i]= 0;
	}
	int max =0;
	for ( i = 0; i < nb_cycles;i++)
	{
		if(liste[i].nb_atomes > max)
			max = liste[i].nb_atomes;
	}
	int min = 0;
	for(i = 1; i < nb_cycles;i++)
	{
		if( position[i] < position[min])
		{
			min = i;
		}
	}
	newposition[0] = min;
	position[min] = max * max;
	for(i = 1; i < nb_cycles;i++)
	{
		//trouver le premier non nul 
			for(j = 0; j < nb_cycles;j++)
			{
				if(position[j] != max *max )
				{
					min = j;
					break;
				}
			}
			for(j = min; j < nb_cycles;j++)
			{
				if( position[j] < position[min])
				{
					min = j;
				}

			}
			newposition[i ] = min;
			position[min] = max * max;

	}

	int **matrice = malloc(nb_cycles * sizeof(int *));
	if( matrice == NULL)
		probleme_memoire();

	for ( i = 0; i < nb_cycles;i++)
	{
		matrice[i] = malloc(m.nb_liaisons  * sizeof(int));
		if( matrice[i] == NULL)
			probleme_memoire();
	}
	int **matrice2 = malloc(nb_cycles * sizeof(int *));
	if( matrice2 == NULL)
		probleme_memoire();

	for ( i = 0; i < nb_cycles;i++)
	{
		matrice2[i] = malloc(m.nb_liaisons  * sizeof(int));
		if( matrice2[i] == NULL)
			probleme_memoire();
	}

	//initialisationd e la matrice 
	for ( i = 0; i < nb_cycles;i++)
	{
		for ( j = 0; j < m.nb_liaisons;j++)
			matrice[i][j] = 0;
	}
	int *chemin;
	int taille,pos ;
	for ( i = 0; i < nb_cycles;i++)
	{
		// je travaille sur le cycle liste[newposition[i]]
		//printf("working on %d \n", liste[newposition[i]].id_cycle);
		taille = liste[newposition[i]].chemin1[0] +  liste[newposition[i]].chemin2[0] + 1;
		chemin = concatener_deux_chemins(liste[newposition[i]]);
		for( j = 0; j < taille ;j++ )
		{
			
			
			pos = position_de_arete(chemin[j],chemin[j + 1],m);
			matrice[i][pos] = 1 ;
		}
		free(chemin);
	}

	//sauvegarde 
	for (i = 0; i < nb_cycles; i++)
	{
		for (j= 0; j < m.nb_liaisons; j++)
		{
			matrice2[i][j] = matrice[i][j];
		}
	}
	//gauss 
	int first,k;
	for (i = 0; i < nb_cycles - 1; i++)
	{
		//trouver le premier un 
		first = -1;
		for (j= 0; j < m.nb_liaisons; j++)
		{
			if(matrice[i][j] == 1)
			{
				first = j;
				break;
			}
		}
		if(first != -1)
		{
			for ( k = i+1; k < nb_cycles; k++)
			{
				if( matrice[k][first] == 1)
				{
					for ( l = 0; l < m.nb_liaisons;l++)
					{
						matrice[k][l] = fonction_xor(matrice[k][l],matrice[i][l]);
					}
				}
			}
		}
	}

	// supprimer tous les cycles independant 
	int somme;
	int base = 0;
	int tab_base[nb_cycles];
	for ( i = 0; i < nb_cycles; i++)
		tab_base[i] = 0;
	for (i = 0; i < nb_cycles; i++)
	{
		somme = 0;
		for (j= 0; j < m.nb_liaisons; j++)
		{
			somme += matrice[i][j];
		}
		if( somme != 0)
		{
			base++;
			tab_base[i] = 1;
		}
			

	}
	// pour tous ceux de la base si la combinaison de deux cycles donne un elimine then le rajouter 

	int *t,maxi;
	for(i = 0; i < nb_cycles-1;i++)
	{
		for(j = i + 1;j < nb_cycles;j++)
		{
			if(tab_base[i] * tab_base[j] == 1) //les deux sont dans la base
			{
				t = produit_xor_matrice(t,matrice2[i],matrice2[j],m.nb_liaisons);
				for(l = 0; l < nb_cycles;l++)
				{
					if(liste[newposition[i]].nb_atomes > liste[newposition[j]].nb_atomes)
						maxi = liste[newposition[i]].nb_atomes;
					else
						maxi = liste[newposition[j]].nb_atomes;
					if(verification_egalite_tableaux(t,matrice2[l],m.nb_liaisons) == 1 && tab_base[l] == 0 && liste[newposition[l]].nb_atomes == maxi)
					{
						tab_base[l] = 1;
						break;
					}
				}
				free(t);
			}
		}
	}

	for(i = 0; i < nb_cycles;i++)
	{
		if(tab_base[i] == 1)
		{
			ajouter_cycle_base(liste[newposition[i]],m);
		}
	}
	//liberation de la memoire de la matrice 
	for(i = 0; i < nb_cycles;i++)
	{
		free(matrice[i]);
		free(matrice2[i]);
	}
		
	free(matrice);
	free(matrice2);


}

void ajouter_cycle_base(cycles c,graphemol m)

{
	if (c.nb_atomes > taille_cycle)
		return;
	if(taille_base == 0)
	{
		labase = malloc((taille_base +1)* sizeof(cycles));
	}
	else
	{
		labase = realloc(labase,(taille_base +1)* sizeof(cycles));
	}
	if( labase == NULL)
			probleme_memoire();
	labase[taille_base] = creer_un_cycle(m,c.sommet,c.arete,c.chemin1,c.chemin2);
	labase[taille_base].id_cycle = taille_base;
	labase[taille_base].pere = c.pere;
	taille_base++;


}	
int *produit_xor_matrice(int *t, int *tab1, int *tab2,int nb)
{
	t= malloc(nb*sizeof(int));
	int i;
	for ( i = 0; i < nb; i++)
		t[i] = fonction_xor(tab1[i],tab2[i]);
	return t;
	
}
//verifie si deux tableaux sont les memes

int verification_egalite_tableaux(int *tab1,int *tab2 , int taille)
{
	
	int i, verif = 1 ;
	for( i = 0; i < taille ; i++)
	{
		if( tab1[i] != tab2[i])
		{
			verif = 0;
			break;
		}
	}
	
	return verif;
	
}

int fonction_xor(int a , int b)
{
	int c;
	if( a == 0)
	{
		if( b ==0)
			c = 0;
		else
			c = 1;

	}
	else
	{
		if( b == 0)
			c = 1;
		else 
			c = 0;
	}
	return c;
}
void affiche_matrice_cycles ( int **matrice , int nb_l,int nb_c)
{
	int i ,j;
	printf("la matrice de cycles \n");
	for (i = 0; i < nb_l; i++)
	{
		for (j= 0; j < nb_c; j++)
		{
			printf("%d ",matrice[i][j] );
		}
		printf("\n");
	}
}
int *concatener_deux_chemins(cycles c)
{
	int s = c.chemin1[0] +  c.chemin2[0] + 1;


	int *chemin  = malloc((s+1) * sizeof(int));
	chemin[0] = c.sommet;
	int i,pos = 1;
	for(i = c.chemin1[0] - 1; i > 0; i--)
	{
		chemin[pos] = c.chemin1[i];
		//printf("%d ", chemin[pos]);
		pos++;
	}
	chemin[pos] = c.arete.A1;
	
	pos++;
	chemin[pos] = c.arete.A2;
	pos++;
	for(i = 0; i < c.chemin2[0] - 1; i++)
	{
		
		chemin[pos] = c.chemin2[i + 1];
		pos++;
	}

	chemin[s] = c.sommet;
	return chemin;
}
int position_de_arete(int sommet1, int sommet2,graphemol m)
{
	int pos = -1;
	int i;

	for ( i = 0; i < m.nb_liaisons; i++)
	{
		if((m.liste_liaisons[i].A1 == sommet1 && m.liste_liaisons[i].A2 == sommet2 )||(m.liste_liaisons[i].A2 == sommet1 && m.liste_liaisons[i].A1 == sommet2 ))
		{
			pos = i;
			break;
		}
	}

	if( pos == -1)
	{
		fprintf(stdout, "impossible de trouver l'arete [ %d - %d] dans le graphe\n", sommet1,sommet2);
		exit(600);
	}
	return pos;
}
cycles *ajouter_un_cycle(cycles *liste, int nb_cycles , cycles c,graphemol m)
{
	//afficher_un_cycle(c,m);
	if( liste == NULL)
	{
		liste = malloc((nb_cycles + 1)* sizeof(cycles));

	}
	else
	{
		liste = realloc(liste, (nb_cycles + 1)* sizeof(cycles));
	}
	if( liste == NULL)
			probleme_memoire();
	liste[nb_cycles].nb_atomes = c.nb_atomes;
	liste[nb_cycles].sommet = c.sommet;
	liste[nb_cycles].pere = c.pere;
	liste[nb_cycles].arete.A1 = c.arete.A1;
	liste[nb_cycles].arete.A2 = c.arete.A2;
	liste[nb_cycles].sommets = malloc(c.nb_atomes * sizeof(int));
	if(liste[nb_cycles].sommets ==  NULL)
		probleme_memoire();
	int i;
	for( i = 0; i < c.nb_atomes; i++)
		liste[nb_cycles].sommets[i] = c.sommets[i];
	liste[nb_cycles].chemin1 = malloc(c.chemin1[0] *sizeof(int));
	if(liste[nb_cycles].chemin1 == NULL)
		probleme_memoire();
	for( i = 0; i < c.chemin1[0]; i++)
		liste[nb_cycles].chemin1[i] = c.chemin1[i];

	liste[nb_cycles].chemin2 = malloc(c.chemin2[0]  *sizeof(int));
	if(liste[nb_cycles].chemin2 == NULL)
		probleme_memoire();
	for( i = 0; i < c.chemin2[0] ; i++)
		liste[nb_cycles].chemin2[i] = c.chemin2[i];
	liste[nb_cycles].id_cycle = nb_cycles + 1;
	return liste;
}
void elimination_feuilles(struct molecule m)
{
	
	graphemol g = conversion_mol_graphe(m);
	int *deja_elimine = malloc(g.nb_atomes * sizeof(int));
	int i;
	for (i = 0; i < g.nb_atomes; i++)
		deja_elimine[i] = 0;
	int *degre = calcul_degre_mol(g);
	g = modification_structure_mol(g,degre,deja_elimine);
	//int res = est_connexe_graphemol(g,deja_elimine);
	int sommets  = g.nb_atomes;
	for ( i = 0;  i < g.nb_atomes;i++)
	{
		//printf("-%d %d\n", i+1,deja_elimine[i]);
		if( deja_elimine[i] == 1)
			sommets--;
	}
	int l;
	if ( sommets > 0 && sommets < 100 )
	{	
		graphemol *liste_connexe = ensemble_connexe_graphemol(g, deja_elimine);
		//printf("Molecule de CHEBI ID :  %d et res =  %d et connexe = %d\n",m.chebi_id ,res,liste_connexe[0].nb_connexe);

		int nb_connexe = liste_connexe[0].nb_connexe;
		int nb_isthmes = 0;
		int *elimine;
		for (i = 0; i < nb_connexe; i++)
		{
			//printf("i = %d\n", i + 1);
			//affiche_graphemol(liste_connexe[i]);
			elimine = malloc(liste_connexe[i].nb_atomes * sizeof(int));
			for ( l = 0;  l < liste_connexe[i].nb_atomes;l++)
			{	
				elimine[l] = 0;
			}

			//printf("magma : %d nb_isthmes = %d\n", i+1 ,nombre_isthmes_graphemol(liste_connexe[i],elimine));
			nb_isthmes += nombre_isthmes_graphemol(liste_connexe[i],elimine);
			free(elimine);
		}

		isthmes *liste_isthmes = malloc(nb_isthmes *sizeof(isthmes));
		if( liste_isthmes == NULL)
			probleme_memoire();

		liste_isthmes = retrouver_tous_isthmes(liste_connexe,liste_isthmes,deja_elimine);

		//for(i = 0; i< nb_isthmes;i++)
		//	printf("%d - %d %d\n",i+1,liste_isthmes[i].l.A1,liste_isthmes[i].l.A2 );
		graphemol h = enlever_tous_isthmes(g,deja_elimine,liste_isthmes,nb_isthmes);

		deja_elimine = realloc(deja_elimine,h.nb_atomes *sizeof(int));
		for (i = 0; i < h.nb_atomes; i++)
			deja_elimine[i] = 0;
		graphemol *t = ensemble_connexe_graphemol(h,deja_elimine);
		//printf("nb total composants bas niveau %d at = %d\n",t[0].nb_connexe,t[0].nb_atomes);	
		//affiche_graphemol(t[0]);
		int p,j;
		//trouver tous les cycles par magma
		for (i = 0; i < t[0].nb_connexe; i++)
		{
			p = t[i].liste_atomes[0];
			for( j = 0; j < nb_connexe;j++)
			{
				if(sommet_dans_basniveau(liste_connexe[j],p))
				{
					//printf("p: %d est dans %d\n",p,j+1 );
					t[i].pere = j;
					break;
				}
			}
			trouver_squelette_cycles(t[i]);
		}
		
		//trouver le composant bas niveauu de chaque cycle de la base
		//printf("taille de la base : %d \n",taille_base );
		for( i = 0; i < taille_base; i++)
		{
			p = labase[i].sommets[0];
			for( j = 0; j < t[0].nb_connexe;j++)
			{
				if(sommet_dans_basniveau(t[j],p))
				{
					//printf("p: %d est dans %d et papy %d\n",p,j+1,t[j].pere );
					labase[i].pere = j;
					break;
				}
			}
		}
		//fichier_base();
		//creation des aretes
		arete_dans_cycle(m);
		arete_dans_cycle_liste(m);
		int dist,dist2;
		ARETE a;
		//affiche_graphemol(h);
		for( i = 0; i < taille_base - 1; i++)
		{
			for( j = i+1; j < taille_base; j++)
			{
				if(labase[i].pere == labase[j].pere) //arete de type 1 ou 2 meme magma
				{
					if(sommets_commun(labase[i],labase[j]) >= 1) // au moins un sommet en commun ---- bleu
					{
						a.id1 = i;
						a.id2 = j;
						a.type = 1;
						a.poids = sommets_commun(labase[i],labase[j]) -1;
						nouvelle_arete(a);
					}
					else
					{// aucune arete en commun
						int r = calcul_LC(labase[i],labase[j],m); //------ rouge
						if( r != -1)
						{	a.id1 = i;
							a.id2 = j;
							a.type = 2;
							a.poids = r;
							nouvelle_arete(a);
							//printf("type 2 arete : ( %d %d %d %d) \n",a.id1+1,a.id2+1,a.type,a.poids );
						}
					}
				}
				else //arete de type 3 magma different 
				{// existe til un chemin de c1 vers c2 
					dist = distance_inter_magma(labase[i],labase[j],m); //----- vert
					if( dist !=-1)
					{//il existe un chemin entre ces deux cycles dans la molecule 
						dist2 = verification_LC(labase[i],labase[j],m,dist);
						if( dist2 != -1)
						{
							a.id1 = i;
							a.id2 = j;
							a.type = 3;
							a.poids = dist;
							nouvelle_arete(a);
						}
					}
				}		
			}
		}
		//liberation de la memoire
		//fichier_dot();
		for (i = 0; i < nb_connexe; i++)
			liberer_graphemol(liste_connexe[i]);
		nb_connexe = t[0].nb_connexe;
		for ( i = 0; i < nb_connexe; i++)
			liberer_graphemol(t[i]);	
		free(liste_connexe);
		free(liste_isthmes);
		liberer_graphemol(h);
		free(t);
		//printf("TAILLE DE LA BASE : %D ",taille_base);
		graphemol sd = conversion_mol_graphe(m);
		if( sd.nb_liaisons > 0)
		{
			for( i = 0; i < sd.nb_liaisons;i++)
				free(arete_liste[i]);
			free(arete_liste);
			free(arete_cycle);
		}
		liberer_graphemol(sd);

		
	}

	
	
	
	free(deja_elimine);
		
	liberer_graphemol(g);
	

	
}
int position_graphemol_arete(graphemol g, int sommet1,int sommet2)
{

	int pos = -1,i;
	for( i = 0; i < g.nb_liaisons;i++)
	{
		if((g.liste_liaisons[i].A1 == sommet1 && g.liste_liaisons[i].A2 == sommet2)||(g.liste_liaisons[i].A1 == sommet2 && g.liste_liaisons[i].A2 == sommet1))
		{
			pos = i;
			break;
		}
	}
	if(pos == -1)
	{
		fprintf(stdout, "liaison non trouve dans la molecule \n" );
		exit(555);
	}
	return pos;
}
void arete_dans_cycle(struct molecule m)
{
	graphemol g = conversion_mol_graphe(m);
	arete_cycle = malloc(g.nb_liaisons * sizeof(int));
	int i,j;
	for ( i = 0; i < g.nb_liaisons;i++ )
		arete_cycle[i] = 0;
	int *chemin,taille;
	for(i = 0; i < taille_base;i++)
	{
		taille = labase[i].chemin1[0] +  labase[i].chemin2[0] + 1;
		
		chemin = concatener_deux_chemins(labase[i]);
		for( j = 0; j < taille ; j++)
		{
			arete_cycle[position_graphemol_arete(g,chemin[j],chemin[j+1])]++;
		}

		free(chemin);
	}

	liberer_graphemol(g);
}

void arete_dans_cycle_liste(struct molecule m)
{
	graphemol g = conversion_mol_graphe(m);
	arete_liste = malloc(g.nb_liaisons *sizeof(int *));
	if(arete_liste == 0)
		probleme_memoire();
	int i;
	for ( i = 0; i < g.nb_liaisons;i++)
	{
		arete_liste[i] = malloc(arete_cycle[i]* sizeof(int));
		if( arete_liste[i] ==  NULL)
			probleme_memoire();
	}
	int *chemin,taille,k,j,pos;
	for ( k = 0; k < g.nb_liaisons;k++)
	{
		pos =0;
		//printf("-- %d :", k);
		for(i = 0; i < taille_base;i++)
		{
			taille = labase[i].chemin1[0] +  labase[i].chemin2[0] + 1;
			
			chemin = concatener_deux_chemins(labase[i]);
			for( j = 0; j < taille ; j++)
			{
				if((g.liste_liaisons[k].A1 == chemin[j] && g.liste_liaisons[k].A2 == chemin[j+1])||(g.liste_liaisons[k].A1 == chemin[j+1] && g.liste_liaisons[k].A2 == chemin[j]))
				{
					arete_liste[position_graphemol_arete(g,chemin[j],chemin[j+1])][pos] = i;
					pos++;
				}
			}

			free(chemin);
		}

	}

	//affichage
	liberer_graphemol(g);

}	
int nombre_de_cycles(int sommet1,int *chemin,int sommet2,graphemol g)
{
	int taille = 0;
	int tableau[taille_base];
	int i,r;
	for( i = 0; i < taille_base;i++)
		tableau[i] = 0;
	if(chemin[0] == 1)
	{
		taille = arete_cycle[position_de_arete(sommet1,sommet2,g)];
	}
	else
	{
		if(arete_cycle[position_de_arete(sommet1,chemin[chemin[0]-1],g)] > 0)
		{
			for( i = 0; i< arete_cycle[position_de_arete(sommet1,chemin[chemin[0]-1],g)] ;i++)
			{
				tableau[arete_liste[position_de_arete(sommet1,chemin[chemin[0]-1],g)][i]] = 1;
			}
		}
		//milieu a faire 
		if(chemin[0]>=3)
		{
			for( r = 1; r <chemin[0]-1;r++)
			{
				if(arete_cycle[position_de_arete(chemin[r],chemin[r+1],g)] > 0)
				{
					for( i = 0; i< arete_cycle[position_de_arete(chemin[r],chemin[r+1],g)] ;i++)
					{
						tableau[arete_liste[position_de_arete(chemin[r],chemin[r+1],g)][i]] = 1;
					}
				}

			}

		}
		if(chemin[0] >= 2)
		{
			if(arete_cycle[position_de_arete(sommet2,chemin[1],g)] > 0)
			{
				for( i = 0; i< arete_cycle[position_de_arete(sommet2,chemin[1],g)] ;i++)
				{
					tableau[arete_liste[position_de_arete(sommet2,chemin[1],g)][i]] = 1;
				}
			}
		}
		for (i = 0; i < taille_base;i++)
			taille += tableau[i];
	}
	return taille;
}
int calcul_LC(cycles a, cycles b,struct molecule m)
{
	
	graphemol g = conversion_mol_graphe(m);
	int res = taille_base * taille_base;
	int i,j;
	int *chemin;
	int valeur,distance = g.nb_liaisons;
	for ( i = 0; i < a.nb_atomes;i++)
	{
		for(j = 0; j < b.nb_atomes;j++)
		{
			chemin = plus_court_chemin(a.sommets[i],b.sommets[j],g);
			
			if(a.sommets[i] != b.sommets[j])
			{
				valeur = nombre_de_cycles(a.sommets[i],chemin,b.sommets[j],g);
				if(valeur <= res && distance >= chemin[0])
				{
					res = valeur;
					distance  = chemin[0];
				}
			}
			free(chemin);
		}
	}
	liberer_graphemol(g);
	if(res < distance_LC )
		res = distance;
	else 
		res  = -1;
	return res;
}


ARETE copier_arete(ARETE a , ARETE b)
{
	a.id1 = b.id1;
	a.id2 = b.id2;
	a.type = b.type;
	a.poids = b.poids;

	return a;
}
GRAPHE_CYCLE construction_graphe_cycles(struct molecule m)
{
	
	elimination_feuilles(m);
	GRAPHE_CYCLE c;
	
	c.nb_sommets = taille_base;
	c.nb_aretes = nb_arete_base;
	int i;
	c.liste_sommets = NULL;
	c.liste_aretes  = NULL;
	c.liste_sommets = malloc( c.nb_sommets * sizeof(SOMMET));
	c.liste_aretes  = malloc(c.nb_aretes * sizeof(ARETE));
	//printf("nb sommets %d et aretes %d\n", c.nb_sommets,c.nb_aretes);
	for( i = 0; i < c.nb_sommets;i++)
	{
		//printf("%d id en cours %d\n",labase[i].id_cycle,labase[i].nb_atomes );
		c.liste_sommets[i].id = labase[i].id_cycle;
		c.liste_sommets[i].poids = labase[i].nb_atomes;

	}
	for( i = 0; i < c.nb_aretes;i++)
	{
		c.liste_aretes[i] = copier_arete(c.liste_aretes[i],base_aretes[i]);
	}
	if(taille_base > 0)
		{
			for( i = 0; i < taille_base;i++)
			{
				liberer_un_cycle(labase[i]);
			}
			free(labase);
			
		}
	if(nb_arete_base > 0)
	{
		free(base_aretes);
	}

	return c;
}
void liberer_graphe_cycles( GRAPHE_CYCLE c)
{
	//int i ; 
	if(c.liste_sommets !=NULL)
		free(c.liste_sommets);
	if(c.liste_aretes !=NULL)
	free(c.liste_aretes);
	
}
int verification_LC(cycles a, cycles b,struct molecule m,int dist)
{
	
	graphemol g = conversion_mol_graphe(m);
	int res = taille_base * taille_base;
	int i,j;
	int *chemin;
	int valeur,distance = g.nb_liaisons;
	for ( i = 0; i < a.nb_atomes;i++)
	{
		for(j = 0; j < b.nb_atomes;j++)
		{
			chemin = plus_court_chemin(a.sommets[i],b.sommets[j],g);
			
			if(a.sommets[i] != b.sommets[j])
			{
				valeur = nombre_de_cycles(a.sommets[i],chemin,b.sommets[j],g);
				if(valeur == 0 && distance >= chemin[0])
				{
					res = valeur;
					distance  = chemin[0];
				}
			}
			free(chemin);
		}
	}
	liberer_graphemol(g);
	if(res == taille_base * taille_base || distance != dist)
		res  = -1;
	return res;
}
void nouvelle_arete(ARETE a)
{
	if(nb_arete_base == 0)
	{
		base_aretes= malloc((nb_arete_base+1)* sizeof(ARETE));
	}
	else
	{
		base_aretes= realloc(base_aretes,(nb_arete_base+1)* sizeof(ARETE));
	}
	if( base_aretes == NULL)
		probleme_memoire();
	base_aretes[nb_arete_base].id1=a.id1;
	base_aretes[nb_arete_base].id2=a.id2;
	base_aretes[nb_arete_base].type=a.type;
	base_aretes[nb_arete_base].poids=a.poids;

	nb_arete_base++;
}
int sommets_commun ( cycles a, cycles b)
{
	int i,j,somme = 0;

	for (i = 0; i < a.nb_atomes; i++)
	{
		for (j = 0; j < b.nb_atomes; j++)
		{
			if(a.sommets[i] == b.sommets[j])
				somme++;
		}
	}


	return somme;
}
int distance_inter_magma(cycles a , cycles b , struct molecule m)
{
	graphemol h = conversion_mol_graphe(m);
	int dist = h.nb_atomes * h.nb_atomes;
	int somme;
	struct liste_voisins * v = construction_voisinage_graphemol(h);
	if(!existe_chaine_graphemol(position_graphemol(h,a.sommet),position_graphemol(h,b.sommet),h,v))
		dist = -1;
	else
	{	
		int i,j,l,*pcc; 
		for( i = 0; i < a.nb_atomes;i++)
		{
			for( j = 0; j < b.nb_atomes;j++)
			{
				pcc = plus_court_chemin(a.sommets[i],b.sommets[j],h);
				if(pcc[0] < dist)
				{
					somme= 0;
					if(pcc[0] > 2)
					{
						somme += arete_cycle[position_graphemol_arete(h,a.sommets[i],pcc[pcc[0]-1])];
						somme += arete_cycle[position_graphemol_arete(h,pcc[1],b.sommets[j])];
						//printf("%d %d \n",a.sommets[i],pcc[pcc[0]-1]);
						//printf("%d %d --- > %d\n",a.sommets[i],pcc[1],arete_cycle[position_graphemol_arete(h,a.sommets[i],pcc[pcc[0]-1])] );
					}
					for( l = 1; l < pcc[0] - 1 ;l++)
						somme += arete_cycle[position_graphemol_arete(h,pcc[l],pcc[l+1])];
					if(somme == 0)
						dist = pcc[0];
				}
					
				free(pcc);
			}
		}
	}
	liberer_memoire_voisins(v,h);
	liberer_graphemol(h);
	return dist;
}
int sommet_dans_basniveau(graphemol t,int p)
{
	int res = 0;
	int i;
	for(i = 0; i < t.nb_atomes;i++)
	{
		if(t.liste_atomes[i] == p)
		{
			res  = 1;
			break;
		}
	}

	return res;
}

void fichier_base()
{
	/*FILE *f;
	char nom[128];
	sprintf(nom,"../CHEBI/%d/%d_%d_%d.base",chebi_id,chebi_id,taille_cycle,distance_LC);
	f = fopen(nom,"w");
	int i;
	fprintf(f, "nb cycles : %d \n",taille_base );
	for ( i = 0; i < taille_base; i++)
		afficher_un_cycle_fichier(labase[i],f);
	fclose(f);*/
}

void afficher_un_cycle_fichier(cycles c,FILE *f)
{
	int i;
	fprintf(f,"id_cycle : %d pere : %d nb atomes %d , sommet = %d et liaison = [%d - %d]\n",c.id_cycle + 1,c.pere+1,c.nb_atomes,c.sommet,c.arete.A1,c.arete.A2 );
	fprintf(f,"Le cycle : %d -",c.sommet);
	for(i  = c.chemin1[0] - 1; i > 0; i--)
		fprintf(f," %d -",c.chemin1[i]);
	fprintf(f," %d  - %d -",c.arete.A1,c.arete.A2);
	for(i  = 0; i < c.chemin2[0] - 1; i++)
		fprintf(f," %d -",c.chemin2[i+1]);
	fprintf(f," %d\n\n",c.sommet);
}


struct molecule construction_matrice_mol(struct molecule m)
{//construction de la matrice de liaison d'une molecule
	
	int i,j;
	if(m.matrice_liaisons == NULL)
	{
		m.matrice_liaisons =  malloc(m.nb_atomes * sizeof(int *));
		
		for(i=0;i< m.nb_atomes;i++) m.matrice_liaisons[i] =  malloc(m.nb_atomes * sizeof(int));
		
		for(i=0;i< m.nb_atomes;i++)
		{
			for(j=0;j< m.nb_atomes;j++)
			 m.matrice_liaisons[i][j] = AUCUNE_LIAISON;
		}

		for(i =0; i< m.nb_liaisons;i++)
		{
			m.matrice_liaisons[m.liste_liaisons[i].A1-1][m.liste_liaisons[i].A2-1]=m.liste_liaisons[i].l_type;
			m.matrice_liaisons[m.liste_liaisons[i].A2-1][m.liste_liaisons[i].A1-1]=m.liste_liaisons[i].l_type;
		}
		

	}

	return m;
	
}


void affiche_matrice(struct molecule m)
{//affcihe la matrice d'une molecule m
	int i,j;
	printf("Affichage de la matrice m \n");
	
	if( m.matrice_liaisons == 	NULL)
	{
		printf("La matrice de cette molecule n'a pas encore eté défini \n");
		return;
	}
	for(i=0;i< m.nb_atomes;i++)
	{
		for(j=0;j< m.nb_atomes;j++)
			printf("%d ",m.matrice_liaisons[i][j] );
		
		printf("\n");
	}
	
}

void affiche_mol(struct molecule M)
{
	int i;
	printf("\nCHEBI ID: %d\n", M.chebi_id);	
	printf("nombre atomes : %d\n", M.nb_atomes);
	printf("nombre de liaisons : %d\n", M.nb_liaisons);
	printf("CHEBI name : ");
		
	for(i=0;i<1024;i++) 
	{
		if ( M.chebi_name[i] !='\n')
			printf("%c", M.chebi_name[i]);
		else
			break;
	}
	printf("\n");
		
	for(i=0;i<M.nb_atomes;i++) printf("%s%d ", atom_name[M.liste_atomes[i]],i + 1);
	printf("\n");
		
	for(i=0;i<M.nb_liaisons;i++) printf("%s%.2d --- %s%.2d %d \n", atom_name[M.liste_atomes[M.liste_liaisons[i].A1 - 1 ]],M.liste_liaisons[i].A1,atom_name[M.liste_atomes[M.liste_liaisons[i].A2 -1]], M.liste_liaisons[i].A2, M.liste_liaisons[i].l_type);
	printf("\n");	
}



//fonctions de similarite

int valeur_absolue(int a)
{
	if ( a > 0)
		return a;
	return -a;
}
int min( int a , int b)
{
	if ( a < b )
		return a;
	return b;
}
int type_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b)
{
	int res = -1;

	int i;
	for ( i = 0; i < c.nb_aretes;i++)
	{
		if((c.liste_aretes[i].id1 == a && c.liste_aretes[i].id2 == b)||(c.liste_aretes[i].id1 == b && c.liste_aretes[i].id2 == a))
		{
			res = c.liste_aretes[i].type;
			break;
		}	
	}
	return res;
}
int poids_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b)
{
	int res = -1;

	int i;
	for ( i = 0; i < c.nb_aretes;i++)
	{
		if((c.liste_aretes[i].id1 == a && c.liste_aretes[i].id2 == b)||(c.liste_aretes[i].id1 == b && c.liste_aretes[i].id2 == a))
		{
			res = c.liste_aretes[i].poids;
			break;
		}	
	}
	return res;
}
GRAPHE_CYCLE construction_graphe_produit(GRAPHE_CYCLE a, GRAPHE_CYCLE b)
{
	GRAPHE_CYCLE c;
	c.liste_aretes = NULL;
	c.liste_sommets = NULL;
	int taille = 0;
	//printf("nb cycles de a %d et nb cycles de b %d\n",a.nb_sommets,b.nb_sommets );
	int i,j;
	for ( i = 0; i < a.nb_sommets; i++)
	{
		//printf("%d :",a.liste_sommets[i].poids );
		for ( j = 0; j < b.nb_sommets; j++)
		{
			if( (float)(valeur_absolue(a.liste_sommets[i].poids - b.liste_sommets[j].poids)) <= 0.2 *min(a.liste_sommets[i].poids ,b.liste_sommets[j].poids))
			{
				c = ajouter_un_sommet(a.liste_sommets[i].id,b.liste_sommets[j].id,taille,c );
				taille++;
				//printf("%d ",b.liste_sommets[j].poids );
			}
		}
	}
	c.nb_sommets = taille;
	int aretes = 0;
	c.liste_aretes = NULL;
	int res1,res2,res3 = 0,res4 = 1;
	for ( i = 0; i < c.nb_sommets - 1; i++)
	{
		for ( j = i + 1; j < c.nb_sommets; j++)
		{
			res3 = 0;
			res4 = 1;
			//printf("le couplet(%d, %d) --- (%d, %d)\n",c.liste_sommets[i].id,c.liste_sommets[i].poids,c.liste_sommets[j].id,c.liste_sommets[j].poids );
			
			if(c.liste_sommets[i].id == c.liste_sommets[j].id)
			{
				res1 = -1;
			}
			else
			{
				res1 = type_arete_graphe_cycle(a,c.liste_sommets[i].id , c.liste_sommets[j].id);
				res3 = poids_arete_graphe_cycle(a,c.liste_sommets[i].id , c.liste_sommets[j].id);
			}	
			if(c.liste_sommets[i].poids == c.liste_sommets[j].poids)
			{
				res2 = -1;
			}
			else
			{
				res2 = type_arete_graphe_cycle(b,c.liste_sommets[i].poids , c.liste_sommets[j].poids);
				res4 = poids_arete_graphe_cycle(b,c.liste_sommets[i].poids , c.liste_sommets[j].poids);
			}	
			if( res1 * res2 >= 0)
			{
				if(res1 == res2 && res3 == res4)
				{
					c = ajouter_une_arete_graphe(i , j,aretes,c);
					aretes++;
					//printf("le couplet(%d, %d) --- (%d, %d) type = %d et poids %d\n",c.liste_sommets[i].id,c.liste_sommets[i].poids,c.liste_sommets[j].id,c.liste_sommets[j].poids,res1,res3 );
			
				}
				/*else 
				{
					if(res1 == res2 && (res1 == 2 || res1 == 3))
					{
						c = ajouter_une_arete_graphe(i , j,aretes,c);
						aretes++;
					}
				}*/
			}
		}
	}
	c.nb_aretes = aretes;
	//printf("taille du graphe produit %d et aretes %d \n", c.nb_sommets, c.nb_aretes);

	return c;
}


GRAPHE_CYCLE ajouter_un_sommet(int a , int b, int taille ,GRAPHE_CYCLE c)
{
	if( c.liste_sommets == NULL)
	{
		c.liste_sommets = malloc((taille + 1)* sizeof(SOMMET));
	}
	else
	{
		c.liste_sommets = realloc(c.liste_sommets,(taille +1)* sizeof(SOMMET));

	}
	if(c.liste_sommets == NULL)
		probleme_memoire();
	c.liste_sommets[taille].id = a;
	c.liste_sommets[taille].poids = b;
	return c;
}

GRAPHE_CYCLE ajouter_une_arete_graphe(int a , int b, int arete ,GRAPHE_CYCLE c)
{
	if( c.liste_aretes == NULL)
	{
		c.liste_aretes = malloc((arete + 1)* sizeof(ARETE));
	}
	else
	{
		c.liste_aretes= realloc(c.liste_aretes,(arete +1)* sizeof(ARETE));

	}
	if(c.liste_aretes == NULL)
		probleme_memoire();
	c.liste_aretes[arete].id1 = a;
	c.liste_aretes[arete].id2 = b;
	c.liste_aretes[arete].poids = 1;
	return c;
}
int ** construction_matrice_produit(GRAPHE_CYCLE p)
{
 int i,j;
 int **matrice = malloc(p.nb_sommets * sizeof(int *));
 if( matrice == NULL)
 	probleme_memoire();
 for(i = 0; i < p.nb_sommets;i++)
 {
 	matrice[i] = malloc(p.nb_sommets * sizeof(int));
 	if(matrice[i] == NULL)
 		probleme_memoire();
 }
 for(i = 0; i < p.nb_sommets;i++)
 {
 	for(j = 0; j < p.nb_sommets;j++)
 	{
 		matrice[i][j] = 0;
 	}
 }
 for(i = 0; i < p.nb_aretes;i++)
 {
 	matrice[p.liste_aretes[i].id1][p.liste_aretes[i].id2] = 1;
 	matrice[p.liste_aretes[i].id2][p.liste_aretes[i].id1] = 1;
 }
 return matrice;
}
void affiche_matrice_produit(int **matrice, int taille)
{
	int i,j;
	for ( i = 0; i < taille ;i++)
	{
 		for(j = 0; j < taille;j++)
 			printf(" %d",matrice[i][j]);
 		printf("\n");
	}
}
void liberation_matrice(int **matrice, int taille)
{
	int i;
	for ( i = 0; i < taille ;i++)
		free(matrice[i]);
	free(matrice);
}

int somme_sommets_aretes( int *clique, GRAPHE_CYCLE produit, int **depart, int sommets)
{

	//printf("i am here \n");
	int somme = 0;
	int i,j;
	
	for( i = 0; i < sommets;i++)
	{
		somme += clique[i];
	}
	for( i = 0; i < sommets - 1 ; i++)
	{
		for( j = i+1; j < sommets; j++)
		{
			
			int id1,id2;
			id1 = produit.liste_sommets[i].id;
			id2 = produit.liste_sommets[j].id;
			if(clique[i] == 1 && clique[j] == 1 && depart[id1][id2] == 1 )
				somme++;
		}
	}
	//printf("%d valeur @@@@@@\n", somme);

	return somme;
}

void calcul_clique(int **matrice,int sommets,int *dans_clique,int taille_clique,int *candidat,int taille_candidat, double date,GRAPHE_CYCLE produit, int **depart)
{	//calcul de la clique max recursif
	

	int i,j;
	
	if( taille_candidat == 0)
	{
		int val1, val2;
		//printf("herer %d. %d \n",taille_clique,taille_clique_max);
		val1 = somme_sommets_aretes(dans_clique_max,produit,depart, sommets);
		val2 = somme_sommets_aretes(dans_clique,produit,depart, sommets);
		if(val2 > val1)
		{
			taille_clique_max = taille_clique;
			for (i = 0 ;  i < sommets ; i++)
				dans_clique_max[i] = dans_clique[i];
		}
		return;
	}
	
	if (taille_candidat + taille_clique <= taille_clique_max)
	{
		//printf("je suis ici \n");
		return;
	}
	
	//else
	//if(taille_candidat == 1)printf("one time\n"); 
	int taille_candidat_temp;
	int *candidat_temp;
	candidat_temp = malloc( sommets * sizeof(int));


	/*if(taille_candidat == 1) {
		for (j = 0 ;  j < sommets ; j ++)
			{	printf("%d ", candidat[j]);}printf(" fin\n");
	}*/
	for (i = 0 ;  i < sommets ; i ++)
	{
		if ( candidat[i] == 1)
		{
			//if(taille_candidat == 1) printf("is %d \n",i);
			candidat[i] = 0;
			dans_clique[i] = 1 ;
			taille_candidat_temp = taille_candidat;
			
			for (j = 0 ;  j < sommets ; j ++)
			{
				candidat_temp[j] = candidat[j]; 
				if ((candidat[j] == 1) && (matrice[i][j] == 0))
				{
					candidat_temp[j] = 0;
					taille_candidat_temp--;	
				}	
			}
			
			taille_candidat_temp--;
			calcul_clique(matrice,sommets,dans_clique,taille_clique + 1,candidat_temp,taille_candidat_temp,date,produit,depart);
			dans_clique[i] = 0;
			candidat[i] = 1;
		}	
		
	}
	free(candidat_temp);
		 
}

	
void la_clique_max( int **matrice, int sommets,int **depart,GRAPHE_CYCLE produit)
{	//Debut calcul de la clique -- Initialisation
	int i;
	int *candidat;
	int *dans_clique;
	
	dans_clique_max = malloc( sommets *sizeof(int));
	if (!dans_clique_max) { fprintf(stderr,"cannot malloc dans_clique_max %d\n",sommets); exit(41); }
	dans_clique = malloc( sommets *sizeof(int));
	if (!dans_clique) { fprintf(stderr,"cannot malloc dans_clique %d\n",sommets); exit(42); }
	candidat = malloc( sommets *sizeof(int));
	if (!candidat) { fprintf(stderr,"cannot malloc candidat %d\n",sommets); exit(43); }
	
	//initialisation 
	for(i = 0; i < sommets ; i++ )
	{
		candidat[i] 	= 1;
		dans_clique[i]	= 0;
		dans_clique_max[i] = 0;
	}
	
	taille_clique_max = 0;

	
	calcul_clique(matrice,sommets,dans_clique,0,candidat,sommets,date_max,produit,depart); // 0 taille de la clique initial  et m.nb_atome = nb sommets candidats
	 
	free(dans_clique);
	free(candidat);
	
}
GRAPHE_CYCLE construction_commun_max(GRAPHE_CYCLE a,GRAPHE_CYCLE p)
{

	GRAPHE_CYCLE c;
	int tab[a.nb_sommets];
	int i, nb_at= 0;
	for(i=0;i < a.nb_sommets ; i++)
	{
		tab[i] = 0;
	}	

	for(i=0;i < p.nb_sommets ; i++)
	{
		
		if(dans_clique_max[i] == 1)
		{
			tab[p.liste_sommets[i].id] = 1;
		}
			
	}
	//free(dans_clique_max);

	for(i=0;i < a.nb_sommets ; i++)
	{
		if(tab[i] ==1)
			nb_at++;
	}
	//printf("na atmes = %d\n",nb_at);
	int j,nb_liaisons = 0;
	for(i=0;i < a.nb_sommets - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j = i+1;j < a.nb_sommets; j++)
			{
				if(tab[j] == 1 && poids_arete_graphe_cycle(a,a.liste_sommets[i].id,a.liste_sommets[j].id) != -1)
						nb_liaisons ++;
			}
		}
	}
	
	c.nb_aretes = nb_liaisons;
	c.nb_sommets = nb_at;
	//printf("na liaisons s = %d\n",nb_liaisons);
	return c;



}

float similarite(GRAPHE_CYCLE a,GRAPHE_CYCLE b)
{

	//int 
	float sim = 0.0;
	if(a.nb_sommets == 0 || b.nb_sommets == 0)
		sim = -1.0;
	else
	{
		last_chrono = chrono();
		GRAPHE_CYCLE produit = construction_graphe_produit(a,b);
		int **matrice_produit = construction_matrice_produit(produit);
		int **depart = construction_matrice_produit(a);
		int i;
		

		la_clique_max(matrice_produit,produit.nb_sommets, depart, produit);
		//affiche_matrice_produit(matrice_produit,produit.nb_sommets);

		for(i = 0; i < a.nb_sommets; i++){
			free(depart[i]);
		}
		free(depart);
		
		//GRAPHE_CYCLE g12 = construction_commun_max(a,produit);
		if(date_max == 0 || (chrono() - last_chrono <= date_max))
		{
		
			GRAPHE_CYCLE g12 = construction_commun_max(a,produit);
			//printf("nb sommets %d et nb_liaison %d\n",g12.nb_sommets,g12.nb_aretes);
	
			float num = (float)((g12.nb_sommets + g12.nb_aretes)*(g12.nb_sommets + g12.nb_aretes));
			float denum = (float)((a.nb_sommets + a.nb_aretes)*(b.nb_aretes + b.nb_sommets));
			sim = num/denum;
		}
		else
		{
			sim = -2;
		}


		free(dans_clique_max);
		taille_clique_max = 0;
		liberation_matrice(matrice_produit,produit.nb_sommets);
		liberer_graphe_cycles(produit);
	}	
	return sim;
}

struct molecule * lecture_fichier_chebi()
{//lecture du fichier chebi.sdf
	
	FILE *F;
	F = fopen("ChEBI_lite.sdf","r");
	
	if ( F == NULL ) 
	{
		 fprintf(stderr,"Cannot open ChEBI_lite.pdf file\n"); 
		 exit(1); 
	}
	init_atom_num();
	int nb_mol, DEB = 0, FIN = NB_MOLECULES;
	struct molecule *M = malloc(NB_MOLECULES*sizeof(struct molecule));
	
	if (M == NULL)
	{
		fprintf(stderr,"Not enough memory for M\n"); 
		exit(3); 
	}
	struct molecule m;
	printf("1. Lecture des molecules : %.3lf s\n",chrono());

	for(nb_mol = DEB ; nb_mol < FIN ; nb_mol++)
	{
		if (nb_mol % 1000 == 0) 
		{ 
			fprintf(stdout,"\r%5d / %d",nb_mol,FIN);
			fflush(stdout); 
		}
		m = lire_molecule_sdf(F);
		M[nb_mol] = m;
	}
	
	fclose(F);
	fprintf(stdout,"\r%5d / %d\n",nb_mol,FIN); 
	printf("Fin de la Lecture des molecules : %.3lf s\n",chrono());
	return M;
	
}
struct molecule lire_molecule_sdf(FILE *F) 
{
	struct molecule M;
	init_atom_num ();
	// On saute l'entête
	ligne_suivante(F);
	ligne_suivante(F);
	ligne_suivante(F);
	// Nombre d'atomes et nombre de liaisons
	M.nb_atomes = lire_entier_3(F);
	M.nb_liaisons = lire_entier_3(F);
//printf("%d %d\n",M.nb_atomes,M.nb_liaisons);
	ligne_suivante(F);
	// Allocation mémoire
	M.liste_atomes   = malloc(M.nb_atomes   * sizeof(int));
	
	if(M.liste_atomes == NULL)
	{
		fprintf(stdout,"Cannot allocate memory  for M.liste_atomes\n");
		exit(55);
	}
	M.liste_liaisons = malloc(M.nb_liaisons * sizeof(struct liaison));
	if(M.liste_liaisons == NULL)
	{
		fprintf(stdout,"Cannot allocate memory  for M.liste_liaisons\n");
		exit(56);
	}
	int l,i,j; 
	M.matrice_liaisons = malloc(M.nb_atomes *sizeof(int *));
	for (l = 0; l < M.nb_atomes; l++)
		M.matrice_liaisons[l] = malloc( M.nb_atomes *sizeof(int));
	M.g.som = NULL;
	M.g.aretes = NULL;
	M.g.matrice_cycles_type = NULL;
	M.g.matrice_cycles_poids = NULL;
	M.nb_hydrogene=0;
	M.g_def = 0;
	M.g.nb_connexe = 0;

	// Lecture du nom des atomes
	int a; for (a=0 ; a<M.nb_atomes ; a++) 
	{
		M.liste_atomes[a] = lire_num_atome(F);
		if(M.liste_atomes[a] == 1) M.nb_hydrogene++;
		ligne_suivante(F);
	}
	
//printf("fin atomes\n");
	// Lecture des liaisons
	for (l=0 ; l<M.nb_liaisons ; l++)
	{
		M.liste_liaisons[l] = lire_liaison(F);
		if(M.liste_liaisons[l].A1 > M.nb_atomes || M.liste_liaisons[l].A2 > M.nb_atomes){ fprintf(stderr,"numero de molecule non valide \n"); exit(4); }
	}
	trouver_la_fin_de_M(F);

	ligne_suivante(F);
	// Lecture du CHEBI ID
	M.chebi_id = lire_chebi_id(F);

	ligne_suivante(F);
	ligne_suivante(F);
	// Lecture du CHEBI NAME
	lire_chebi_name(F,&M);
	// Lit la fin de la molecule ($$$$)
	lire_fin_molecule(F);
	//remplissage de la matrice liaison
	
	for (i=0 ; i<M.nb_atomes ; i++)
	{
		for (j=0 ; j<M.nb_atomes; j++)
			M.matrice_liaisons[i][j] = 0;
	}
	for (l=0 ; l<M.nb_liaisons ; l++)
	{
		//printf("%d %d\n",M.liste_liaisons[l].A1,M.liste_liaisons[l].A2 );
		M.matrice_liaisons[M.liste_liaisons[l].A1 -1][M.liste_liaisons[l].A2 -1] = M.liste_liaisons[l].l_type;
		M.matrice_liaisons[M.liste_liaisons[l].A2 -1][M.liste_liaisons[l].A1 -1] = M.liste_liaisons[l].l_type;
	}

	return M;
}
int position_M( int g1_chebi,struct molecule *M)
{ // trouve la poition d'une molecule de chebi g1_chebi dans M
	int i;
	
	for(i=0;i< NB_MOLECULES ;i++)
	{
		if(M[i].chebi_id == g1_chebi)
			return i;
	}
	if(i == NB_MOLECULES)
	{ 
		fprintf(stderr,"numero de chebi non present dans la base \n");
		exit(5); 
	}
	
	return 0;
}
