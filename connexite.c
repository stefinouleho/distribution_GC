#include "modelisation.h"


struct liste_voisins*  construction_voisinage_graphemol(graphemol M){
	
	struct liste_voisins *voisins;
	int i,j,k,a1, a2;
	int nb_sommets = M.nb_atomes;
	int tableau_nb_voisins[nb_sommets];
	//initialisation du tableau à 0
	for(i = 0; i < nb_sommets; i++)
		tableau_nb_voisins[i] = 0;
	//remplissage
	int pos1,pos2;	
	for(i = 0; i < M.nb_liaisons; i++)
	{
		
		a1 = M.liste_liaisons[i].A1;
		a2 = M.liste_liaisons[i].A2;
		pos1= position_graphemol(M,a1);
		pos2= position_graphemol(M,a2);
		tableau_nb_voisins[ pos1]++;
		tableau_nb_voisins[ pos2]++;
		
	
	}
	//allocation memoire
	voisins = malloc( nb_sommets * sizeof(struct liste_voisins));
	if( voisins == NULL)
		probleme_memoire();
	
	for(i = 0; i < nb_sommets; i++)
	{
		//printf("i = %d; atome = %d et voisins = %d \n",i,M.liste_atomes[i], tableau_nb_voisins[i]);
		voisins[i].id_atome = M.liste_atomes[i] ;
		voisins[i].nb_voisins = tableau_nb_voisins[i];
		if( voisins[i].nb_voisins  > 0 )
		{
		
			voisins[i].id_voisins = malloc(voisins[i].nb_voisins * sizeof(int));
			
			if(voisins[i].id_voisins == NULL)
				probleme_memoire();
			k = 0;
		
			for(j = 0; j < M.nb_liaisons; j++)
			{
				if( M.liste_liaisons[j].A1 == voisins[i].id_atome )
				{
					voisins[i].id_voisins[k] = M.liste_liaisons[j].A2;
					//printf("%d %d\n", k+1,voisins[i].id_voisins[k]);
					k++;
				}
				
				else if ( M.liste_liaisons[j].A2 == voisins[i].id_atome )
				{
					voisins[i].id_voisins[k] = M.liste_liaisons[j].A1;
					//printf("%d %d\n", k+1,voisins[i].id_voisins[k]);
					k++;
				}
				
				
			}
		} 

	}
	//affichage_liste_voisinage_graphemol(voisins,M);
	return voisins;
}


//affichage liste voisinage de la molecule

void affichage_liste_voisinage_molecule(struct liste_voisins* voisins,struct molecule M)
{
	printf("Affichage du voisinage de la molecule \nNombre d'atomes est de %d et nb liaisons : %d\n",M.nb_atomes,M.nb_liaisons);
	
	int i,j; 
	for( i = 0; i < M.nb_atomes; i++)
	{
			printf("%s%.2d ---> ",atom_name[M.liste_atomes[i]],i+1);
			
			if(voisins[i].nb_voisins > 0)
			{
				for(j = 0; j < voisins[i].nb_voisins; j++)
				{
					printf("%s%.2d ",atom_name[M.liste_atomes[voisins[i].id_voisins[j]-1]],voisins[i].id_voisins[j]);
				}
			}
			else
			{
				printf(" - ");
			}
			printf(" \n");
	}
	
}

void affichage_liste_voisinage_graphemol(struct liste_voisins* voisins,graphemol M)
{
	printf("Affichage du voisinage de la molecule \nNombre d'atomes est de %d et nb liaisons : %d\n",M.nb_atomes,M.nb_liaisons);
	
	int i,j,pos,k; 
	for( i = 0; i < M.nb_atomes; i++)
	{
			printf("%s%.2d ---> ",atom_name[M.type_atomes[i]],M.liste_atomes[i]);
			pos=-1;
			if(voisins[i].nb_voisins > 0)
			{
				for(j = 0; j < voisins[i].nb_voisins; j++)
				{
					//printf("%s%.2d ",atom_name[M.type_atomes[M.type_atomes[voisins[i].id_voisins[j]]]],voisins[i].id_voisins[j]);
					for(k = 0; k < M.nb_atomes;k++)
					{
						if(M.liste_atomes[k] == voisins[i].id_voisins[j])
						{
							pos = k;
							break;
						}
					}
					if(pos == -1)
					{
						fprintf(stdout, "Impossible de trouver la position de cet atome dans la molecule\n");
						exit(457);
					}
					printf("%s%.2d ",atom_name[M.type_atomes[pos]],voisins[i].id_voisins[j]);
				}
			}
			else
			{
				printf(" - ");
			}
			printf(" \n");
	}
	
}


int position_graphemol(graphemol g,int element)
{
	int i ;
	int pos = -1;
	for ( i = 0; i < g.nb_atomes;i++)
	{
		if(g.liste_atomes[i]== element)
		{	pos = i;
			break;
		}
	}
	if( pos == -1)
	{
		fprintf(stdout, "Impossible de trouver la position de cet atome dans la molecule\n");
		exit(458);
	}
	return pos;
}
//verifie si dans le graphe de la molecule m il existe une chaine entre les sommets en position i et j.( 1 si oui et 0 sinon)

int existe_chaine_graphemol(int i, int j , graphemol m,struct liste_voisins * v)
{
	int resultat = 0,k,l,somme,somme_new,nb_valeurs;
	int trouve[m.nb_atomes];
	int pos;
	//initialise à distance non atteignable pour tous les autres atomes
	for (k = 0; k < m.nb_atomes; k++)
		trouve[k] = -1;
	//metre à jour les voisins les plus proches ( distance 1)
	for (k = 0; k < v[i].nb_voisins; k++)
	{
		pos = position_graphemol(m,v[i].id_voisins[k]);
		trouve[pos] = 1;
	}
		
	trouve[i] = 0;
	//parcours en largeur
	somme = (v[i].nb_voisins - m.nb_atomes);  //nombre de d'atomes pas encore atteint
	somme_new = 0;
	int iteration = 1;
	while ( trouve[j] == -1 && somme != somme_new )
	{
		somme_new = somme;
		nb_valeurs = 0;
		for (k = 0; k < m.nb_atomes; k++)
		{
			if(trouve[k] == iteration)
			{
				//printf("iteration = %d sommet = %d\n",iteration ,k +1);
				for (l = 0; l < v[k].nb_voisins; l++)
				{
					pos = position_graphemol(m,v[k].id_voisins[l]);
					if( trouve[pos] == -1)
					{
						trouve[pos ] = trouve[k] + 1;
						nb_valeurs++;
					}
						
				}
			}
		}
		somme = somme_new + nb_valeurs;
		iteration = iteration + 1;
	}
	if( trouve[j] != -1)
		resultat = 1;

	return resultat;
}
//verifie si le graphe est connexe, un booleen indiquant si m est connexe ou non ( 1 si oui et 0 sinon)


int est_connexe_graphemol2(graphemol m, int *deja_elimine)
{
	//affiche_mol(m);
	int resultat  = 1,i,j,chaine;
	//construction de la liste des voisins
	struct liste_voisins *v =  construction_voisinage_graphemol(m);

	for (i = 0; i < m.nb_atomes - 1; i++)
	{
		//trouver s'il existe un chemin de l'atome en position i vers l'atome en position jj
		for (j = i + 1; j < m.nb_atomes; j++)
		{
			//verifie s'il existe une chaine entre i et j dans m
			if(deja_elimine[i] == 0 && deja_elimine[j] == 0)
			{
				chaine = existe_chaine_graphemol(i,j,m,v);
				//printf("existe chaine %d %d %d\n", m.liste_atomes[i],m.liste_atomes[j],chaine);
				if(!chaine)
				{
					resultat = 0;
					break;
				}

			}
			
		}
		if(j != m.nb_atomes && !chaine)
			break;
	}
	for (i = 0; i < m.nb_atomes; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);

	return resultat;
}


int est_connexe_graphemol(graphemol m, int *deja_elimine)
{
	//affiche_mol(m);
	int resultat  = 1,i,j,chaine;
	//construction de la liste des voisins
	struct liste_voisins *v =  construction_voisinage_graphemol(m);

	for (i = 0; i < m.nb_atomes - 1; i++)
	{
		//trouver s'il existe un chemin de l'atome en position i vers l'atome en position jj
		for (j = i + 1; j < m.nb_atomes; j++)
		{
			//verifie s'il existe une chaine entre i et j dans m
			//printf(" c %d %d %d %d \n", m.liste_atomes[i],m.liste_atomes[j],deja_elimine[i],deja_elimine[j]);
			if(deja_elimine[i] == 0 && deja_elimine[j] == 0)
			{
				chaine = existe_chaine_graphemol(i,j,m,v);
				//printf("existe chaine %d %d %d\n", m.liste_atomes[i],m.liste_atomes[j],chaine);
				if(!chaine)
				{
					resultat = 0;
					break;
				}

			}
			
		}
		if(j != m.nb_atomes && !chaine)
			break;
	}
	for (i = 0; i < m.nb_atomes; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);

	return resultat;
}
//trouve tous les sommets atteignables a partir du sommet i

int * sommets_atteignable_graphemol(int i,graphemol m,struct liste_voisins * v,int *deja_elimine)
{
	int k,l,somme,somme_new,nb_valeurs,pos;
	int *trouve = malloc(m.nb_atomes* sizeof(int));

	if( trouve == NULL)
		probleme_memoire();
	
	//initialise à distance non atteignable pour tous les autres atomes
	for (k = 0; k < m.nb_atomes; k++)
		trouve[k] = -1;
	//metre à jour les voisins les plus proches ( distance 1)
	for (k = 0; k < v[i].nb_voisins; k++)
	{
		pos = position_graphemol(m,v[i].id_voisins[k]);
		trouve[pos] = 1;
	}
		
	trouve[i] = 0;
	//parcours en largeur
	int total = 0;
	for ( k = 0; k < m.nb_atomes;k++)
	{
		if(deja_elimine[i] == 1)
			total++;
	}
	somme = (v[i].nb_voisins - m.nb_atomes + total);  //nombre de d'atomes pas encore atteint
	somme_new = 0;
	int iteration = 1;
	while (somme != somme_new )
	{
		somme_new = somme;
		nb_valeurs = 0;
		for (k = 0; k < m.nb_atomes; k++)
		{
			if(trouve[k] == iteration)
			{
				//printf("iteration = %d sommet = %d\n",iteration ,k +1);
				for (l = 0; l < v[k].nb_voisins; l++)
				{
					pos = position_graphemol(m,v[k].id_voisins[l]);
					if( trouve[pos ] == -1)
					{
						trouve[pos] = trouve[k] + 1;
						nb_valeurs ++;
					}
						
				}
			}
		}
		somme = somme_new + nb_valeurs;
		iteration = iteration + 1;
	}
	return trouve;

}
//liberation liste_voisins


void liberation_liste_voisins_graphemol(struct liste_voisins *v,graphemol m)
{
	int i;
//liberation memoire liste des voisins 
	for (i = 0; i < m.nb_atomes; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);
}
void liberation_graphe(struct graphe g)
{
	int i;
	if( g.som != NULL)
		free(g.som);
	if( g.aretes != NULL)
		free(g.aretes);
	if (g.matrice_cycles_type != NULL)
	{
		for (i = 0; i < g.nb_sommets; i++)
			free(g.matrice_cycles_type[i]);
		free(g.matrice_cycles_type);
	}
	if (g.matrice_cycles_poids != NULL)
	{
		for (i = 0; i < g.nb_sommets; i++)
			free(g.matrice_cycles_poids[i]);
		free(g.matrice_cycles_poids);
	}

}


void probleme_memoire()
{
	fprintf(stdout,"Cannot allocate memory\n");
		exit(145);
}
//construit tous les sous graphe induits connexes de la molecule


graphemol * ensemble_connexe_graphemol(graphemol m , int *deja_elimine)
{
	graphemol *g = NULL;
	int *degre = calcul_degre_mol(m);
	int marque[m.nb_atomes];
	int i,j,position,somme = m.nb_atomes;
	for ( i = 0; i < m.nb_atomes; i++)
	{
		if( deja_elimine[i] == 1)
			somme--;
	}

	//affiche_mol(m);
	int nb_connexe = 0;
	struct liste_voisins *v =  construction_voisinage_graphemol(m);
	for ( i = 0; i < m.nb_atomes; i++)
	{
		if(deja_elimine[i] == 0)
			marque[i] = 1;
		else
			marque[i] = 0;
	}
	
	int *t;
	int sommets;
	while( somme > 0)
	{
		nb_connexe++;
		for ( i = 0; i < m.nb_atomes; i++)
		{
			if(marque[i] == 1)
			{
				position = i;
				break;
			}
		}
		//printf("position = %d \n", position);
		t  = sommets_atteignable_graphemol(position,m,v,deja_elimine);

		sommets = 0;
		for ( i = 0; i < m.nb_atomes; i++)
		{	
			//printf("%d ",t[i] );
			if(t[i] >= 0 && marque[i] == 1)
			{
				marque[i] = 0;
				sommets = sommets + 1;
			}

		}
		
		if (g == NULL && nb_connexe == 1)
		{
			//printf("nb connexes %d\n",nb_connexe );
			g = malloc(nb_connexe * sizeof(graphemol));
			if( g == NULL)
				probleme_memoire();


		
		}
		else
		{	
			//printf("nb connexes %d\n",nb_connexe );
			g = realloc(g,nb_connexe * sizeof(graphemol));
			if( g == NULL)
				probleme_memoire();
		}
		
		//mise a jour de la composante 
		g[nb_connexe - 1].nb_atomes = sommets;
		g[nb_connexe - 1].nb_liaisons = 0;
		g[nb_connexe - 1].pere = m.pere + 1;
		int aretes = 0;

		//compter les aretes
		int pos1,pos2;
		for ( i = 0; i < m.nb_liaisons;i++)
		{
			pos1 = position_graphemol(m,m.liste_liaisons[i].A1);
			pos2 = position_graphemol(m,m.liste_liaisons[i].A2);
			if(t[pos1] >= 0 || t[pos2] >= 0)
				aretes++;
		}
		//stockage des sommets
		//printf("connexe , %d\n",nb_connexe );
		g[nb_connexe - 1].liste_atomes = malloc(sommets *sizeof(int));
		g[nb_connexe - 1].type_atomes = malloc(sommets *sizeof(int));

		if(g[nb_connexe - 1].liste_atomes== NULL)
			probleme_memoire();
		if(g[nb_connexe - 1].type_atomes== NULL)
			probleme_memoire();
		i = 0;
		int compteur = 0;
		while(i < m.nb_atomes)
		{
			
			if(t[i] >=0)
			{
				
				g[nb_connexe - 1].liste_atomes[compteur] = m.liste_atomes[i];
				g[nb_connexe - 1].type_atomes[compteur] = m.type_atomes[i];
				compteur++;
			}
			
			i = i + 1 ;
		}
		//stockage des aretes
		g[nb_connexe - 1].nb_liaisons = aretes; 
		g[nb_connexe - 1].liste_liaisons = malloc(aretes *sizeof(struct liaison));
		if(g[nb_connexe - 1].liste_liaisons == NULL)
			probleme_memoire();
		//stockage du type des liaisons
		g[nb_connexe - 1].nb_connexe = nb_connexe;
		g[nb_connexe - 1].matrice_liaisons = malloc(g[nb_connexe - 1].nb_atomes *sizeof(int *));
		for ( i = 0; i < g[nb_connexe - 1].nb_atomes; i++)
		{
			g[nb_connexe - 1].matrice_liaisons[i] = malloc(g[nb_connexe - 1].nb_atomes * sizeof(int));
			if(g[nb_connexe - 1].matrice_liaisons[i] == NULL)
				probleme_memoire();	
		}
		

		for ( i = 0; i < g[nb_connexe - 1].nb_atomes; i++)
		{
			for ( j = 0; j < g[nb_connexe - 1].nb_atomes; j++)
				g[nb_connexe - 1].matrice_liaisons[i][j] = 0;
		}


		i = 0;
		compteur = 0;
		while ( i < m.nb_liaisons)
		{
			pos1 = position_graphemol(m,m.liste_liaisons[i].A1);
			pos2 = position_graphemol(m,m.liste_liaisons[i].A2);
			if(t[pos1] >= 0 || t[pos2] >= 0)
			{
				g[nb_connexe - 1].liste_liaisons[compteur].A1 = m.liste_liaisons[i].A1;
				g[nb_connexe - 1].liste_liaisons[compteur].A2 = m.liste_liaisons[i].A2;
				g[nb_connexe - 1].liste_liaisons[compteur].l_type = m.liste_liaisons[i].l_type;
				compteur++;

			}
			i = i + 1 ;
		}
		for ( i = 0; i < g[nb_connexe -1].nb_liaisons; i++)
		{
			pos1 = position_graphemol( g[nb_connexe - 1],g[nb_connexe - 1].liste_liaisons[i].A1);
			pos2 = position_graphemol( g[nb_connexe - 1],g[nb_connexe - 1].liste_liaisons[i].A2);

			g[nb_connexe - 1].matrice_liaisons[pos1][pos2]=g[nb_connexe - 1].liste_liaisons[i].l_type;
			g[nb_connexe - 1].matrice_liaisons[pos2][pos1]=g[nb_connexe - 1].liste_liaisons[i].l_type;
		}

		



		//recalcul de la somme
		somme = 0;
		for ( i = 0; i < m.nb_atomes; i++)
			somme += marque[i];
		//printf("somme = %d\n", somme);

		free(t);
		
	}
	free(degre);
	liberation_liste_voisins_graphemol(v,m);

	g[0].nb_connexe = nb_connexe ;
	return g;
}
graphemol enlever_une_arete(graphemol g, struct liaison l)
{
	int i,pos1,pos2;
	pos1 = position_graphemol(g,l.A1);
	pos2 = position_graphemol(g,l.A2);

	g.matrice_liaisons[pos1][pos2]= 0;
	g.matrice_liaisons[pos2][pos1]= 0;

	//stockage
	struct liaison store[g.nb_liaisons];
	for ( i = 0; i < g.nb_liaisons; i++)
	{
		store[i].A1 = g.liste_liaisons[i].A1;
		store[i].A2 = g.liste_liaisons[i].A2;
		store[i].l_type = g.liste_liaisons[i].l_type;
	}

	g.nb_liaisons = g.nb_liaisons - 1;

	int compteur = 0; 
	for ( i = 0; i < g.nb_liaisons + 1; i++)
	{	
		if(store[i].A1 != l.A1 || store[i].A2 != l.A2 )
		{
			g.liste_liaisons[compteur].A1 = store[i].A1;
			g.liste_liaisons[compteur].A2 = store[i].A2;
			g.liste_liaisons[compteur].l_type = store[i].l_type;
			compteur++;
		}
	}
	g.liste_liaisons[g.nb_liaisons ].A1 =0;
	g.liste_liaisons[g.nb_liaisons ].A2 =0;
	g.liste_liaisons[g.nb_liaisons ].l_type =0;

	return g;
}
graphemol ajouter_une_arete(graphemol g, struct liaison l)
{
	int i,p,pos1=-1,pos2=-1;
	for( p = 0; p < g.nb_atomes;p++)
	{
		if(g.liste_atomes[p] == l.A1)
		{
			pos1=p;
			break;
		}
	}

	for( p = 0;p < g.nb_atomes;p++)
	{
		if(g.liste_atomes[p] == l.A2)
		{
			pos2=p;
			break;
		}
	}
	if(pos1 == -1 || pos2 ==-1)
	{
		fprintf(stdout, "La position de l'atome dans la molécule n'a pas été trouvé\n");
		exit(234);
	}	

	g.matrice_liaisons[pos1][pos2]= l.l_type;
	g.matrice_liaisons[pos2][pos1]= l.l_type;
	//stockage
	g.nb_liaisons = g.nb_liaisons + 1;
	for ( i = 0; i < g.nb_liaisons ; i++)
	{
		//printf("--%d %d %d ",i, g.liste_liaisons[i].A1,g.liste_liaisons[i].A2);
		if(g.liste_liaisons[i].A1 == 0 && g.liste_liaisons[i].A2 == 0)
		{
			g.liste_liaisons[i ].A1 =l.A1;
			g.liste_liaisons[i].A2 =l.A2;
			g.liste_liaisons[i ].l_type =l.l_type;
		}
	}
	
	

	return g;
}
int retirer_liaison_connexe(graphemol g, struct liaison l, int *deja_elimine)//on retire la liaison l et on verifie que le graphe reste connexe
{//on retourne 0 si le graphe reste connexe et 1 sinon
	
	int resultat = 0;


	g = enlever_une_arete(g,l);
	
	resultat = 1 - est_connexe_graphemol(g,deja_elimine);
	g = ajouter_une_arete(g,l);
	
	

	return resultat;
}

int nombre_isthmes_graphemol(graphemol g, int *deja_elimine)
{
	int nb_isthmes = 0,i;
	i = 0; 

	while ( i  < g.nb_liaisons)
	{
		i++;
		//printf("%d %d %d \n",g.liste_liaisons[0].A1,g.liste_liaisons[0].A2,nb_isthmes);
		nb_isthmes += retirer_liaison_connexe(g,g.liste_liaisons[0],deja_elimine);
	}


	return nb_isthmes;
}


isthmes *retrouver_tous_isthmes(graphemol *liste_connexe,isthmes *liste, int *deja_elimine)
{
	int nb_connexe = liste_connexe[0].nb_connexe;
	isthmes p;
	int nb_isthmes = 0,i,j,l;
	int *elimine;
	for ( j = 0;  j < nb_connexe ;  j++)
	{
		i = 0; 
		//printf("j = %d nb liaisons = %d\n",j, liste_connexe[j].nb_liaisons );

		elimine = malloc(liste_connexe[j].nb_atomes * sizeof(int));
		for ( l = 0;  l < liste_connexe[j].nb_atomes;l++)
		{	
				elimine[l] = 0;
		}
		while ( i  < liste_connexe[j].nb_liaisons)
		{
			p.l.A1 = liste_connexe[j].liste_liaisons[0].A1;
			p.l.A2 = liste_connexe[j].liste_liaisons[0].A2;
			p.l.l_type = liste_connexe[j].liste_liaisons[0].l_type;
			p.id_composant = liste_connexe[j].pere;
			//printf("j: %d %d %d\n",liste_connexe[j].liste_liaisons[0].A1,liste_connexe[j].liste_liaisons[0].A2, retirer_liaison_connexe(liste_connexe[j],liste_connexe[j].liste_liaisons[0],deja_elimine));
			if(retirer_liaison_connexe(liste_connexe[j],liste_connexe[j].liste_liaisons[0],elimine) == 1)
			{
				liste[nb_isthmes].l.A1 = p.l.A1;
				liste[nb_isthmes].l.A2 = p.l.A2;
				liste[nb_isthmes].l.l_type = p.l.l_type;
				liste[nb_isthmes].id_composant = p.id_composant;
				//printf("i: %d %d \n",liste[nb_isthmes].l.A1,liste[nb_isthmes].l.A2 );  
				nb_isthmes++;
			}
			i++;
		}
		free(elimine);
	}
	//affiche_liste_isthmes(liste,nb_isthmes);
	return liste;
}
void affiche_liste_isthmes ( isthmes *l, int nb)
{
	int i ; 
	printf("Affichage de la liste des isthmes de la molécule \n");
	for (i = 0; i < nb; i++)
	{
		printf("%d [%.2d - %.2d] type = %d et id_composant = %d\n",i ,l[i].l.A1,l[i].l.A2,l[i].l.l_type,l[i].id_composant );
	}
}

graphemol enlever_tous_isthmes(graphemol g,int *deja_elimine,isthmes *liste_isthmes,int nb_isthmes)
{

	graphemol h;
	int *degre = calcul_degre_mol(g);
	int atomes = 0,i,j,pos,pos1;
	for (i = 0; i < nb_isthmes; i++)
	{
		//printf("%d %d ---- %d %d\n", position_graphemol(g,liste_isthmes[i].l.A1),degre[position_graphemol(g,liste_isthmes[i].l.A1)],position_graphemol(g,liste_isthmes[i].l.A2),degre[position_graphemol(g,liste_isthmes[i].l.A2)]);
		pos = position_graphemol(g,liste_isthmes[i].l.A1);
		degre[pos]--;
		pos = position_graphemol(g,liste_isthmes[i].l.A2);
		degre[pos]--;
	}
	for(i = 0; i < g.nb_atomes; i++)
	{
		if(degre[i] > 0)
			atomes++;
	}
	h.nb_atomes = atomes;
	h.nb_liaisons = g.nb_liaisons - nb_isthmes;
	h.liste_atomes = malloc(h.nb_atomes * sizeof(int));
	if(h.liste_atomes == NULL)
		probleme_memoire();
	h.type_atomes = malloc(h.nb_atomes * sizeof(int));
	if(h.type_atomes == NULL)
		probleme_memoire();
	i = 0;

	while ( i < h.nb_atomes)
	{
		for(j = 0; j <g.nb_atomes;j++ )
		{
			if(degre[j] > 0)
			{
				h.liste_atomes[i] = g.liste_atomes[j];
				h.type_atomes[i] = g.type_atomes[j];
				i++;
			}
		}
	}
	h.liste_liaisons = malloc(h.nb_liaisons *sizeof(struct liaison));
	int compteur = 0 ;
	for( i  = 0; i < nb_isthmes; i++)
	{
		for(j = 0;  j < g.nb_liaisons;j++)
		{
			if((g.liste_liaisons[j].A1 == liste_isthmes[i].l.A1 )&&(g.liste_liaisons[j].A2 == liste_isthmes[i].l.A2))
			{
				g.liste_liaisons[j].A1 = 0;
				g.liste_liaisons[j].A2 = 0;
				break;
			}
		}
	}
	for(j = 0;  j < g.nb_liaisons;j++)
	{
		if((g.liste_liaisons[j].A1 != 0 )&&(g.liste_liaisons[j].A2 != 0))
		{
			h.liste_liaisons[compteur].A1 = g.liste_liaisons[j].A1 ;
			h.liste_liaisons[compteur].A2 = g.liste_liaisons[j].A2;
			h.liste_liaisons[compteur].l_type = g.liste_liaisons[j].l_type;
			compteur++;
		}
	}
	h.matrice_liaisons = malloc(h.nb_atomes * sizeof(int *));
	
	for ( i = 0; i < h.nb_atomes; i++)
		h.matrice_liaisons[i] = malloc(h.nb_atomes  *sizeof(int));
	for ( i = 0; i < h.nb_atomes; i++)
	{
		for ( j = 0; j < h.nb_atomes; j++)
		{
			h.matrice_liaisons[i][j] = 0;
		}
	}

	for(j = 0;  j < h.nb_liaisons;j++)
	{
		pos = position_graphemol(h,h.liste_liaisons[j].A1);
		pos1 = position_graphemol(h,h.liste_liaisons[j].A2);
		h.matrice_liaisons[pos][pos1] = h.liste_liaisons[j].l_type;
		h.matrice_liaisons[pos1][pos] = h.liste_liaisons[j].l_type;
	}
	h.nb_connexe  = g.nb_connexe;
	h.pere  = g.pere + 1;
	free(degre);
	return h;
}

