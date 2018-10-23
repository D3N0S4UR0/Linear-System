#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"

int main()
{
    int qtdLinhas, i;
    double **mat;
    FILE *arq;

    res = (Lista*)malloc(sizeof(Lista));
    vars = (Lista*)malloc(sizeof(Lista));

    res->inicio = NULL;
    vars->inicio = NULL;

    arq = fopen("infos.txt", "r");
    if(!arq)
        return 1;

    fscanf(arq, "%d\n", &qtdLinhas);

    mat = inicializaMatriz(qtdLinhas);

    lerDados(mat, qtdLinhas, arq);
    resolve(mat, qtdLinhas);
}

/**
    Lê os dados do arquivo e os armazena nos devidos lugares
    @param mat: a matriz onde os dados devem ser armazenados
    @param n: a ordem da matriz onde os dados devem ser armazenados
    @param file_pointer: o ponteiro do arquivo a ser lido
*/
void lerDados(double **mat, int n, FILE *file_pointer)
{
    double *d;
    long bSeek = 0;
    char *c;
    char f;
    int i, k, t;
    c = (char *)malloc(sizeof(char));
    d = (double *) malloc(sizeof(double));
    for(i=0; i<n; i++)
    {
        *c='@';
        while(*c != '=' && *c != '\n')
        {
            if(*c != '@')
                fseek(file_pointer, -sizeof(char), SEEK_CUR);

            *d = 0;
            fscanf(file_pointer, "%lf", d);
            if(*d == 0)
            {
                fseek(file_pointer, -sizeof(char), SEEK_CUR);
                fscanf(file_pointer, "%c", c);
                if(*c!=85 && *c!=10)
                {
                    if(*c=='-')
                    *d = -1;
                    else
                    *d = 1;
                }
                else
                    *d=1;
            }
            fscanf(file_pointer, "%c", c);
            if(*c =='=')
                break;
            k = tem(vars, (void*) c) ;
            if(k != 0)
            {
                *(*(mat + i) + k-1) = *d;
            }
            else
            {
                k = insira(vars, (void*) c, sizeof(char));
                *(*(mat + i) + k) = *d;
            }
            fscanf(file_pointer, "%c", c);
        }
            fscanf(file_pointer, "%lf", d);
            insira(res, (void*)d, sizeof(double));
    }
    free(c);
    free(d);
}

/**
    verifica se há algum dado na lista coincidente com o valor guardado no endereço apontado pelo ponteiro *inf,
    caso haja retorna sua posição para facilitar a busca do dado
    @param lis: o endereço correspondente a lista onde se deve buscar a informação
    @param inf: o endereço genérico correspondente à informação a ser buscada
    @return a posição em que encontrou o valor na lista
*/
int tem(Lista* lis, void* inf)
{
    int i = 0;
    char *c, *s;

    c = (char*)malloc(sizeof(char));
    s = (char*)malloc(sizeof(char));

    c = (char*) inf;
    No* atual = lis->inicio;
    while(atual != NULL)
    {
        s = (char*) atual->info;
        i++;
        if(*s == *c)
            return i;

        atual = atual->prox;
    }
    return 0;
}

/**
    Já com todos os dados em seus devidos lugares, o método resolve o sistema utilizando o método de Cramer
    @param mat: a matriz dos coeficientes do sistema
    @param n: a ordem da matriz mat
    @return
*/
void resolve(double** mat, int n)
{
    double D, d1, Dx;
    double ** aux;
    int i, j, k;
    D = determinante(mat, n);


    aux =  inicializaMatriz(n);

    if(D == 0)
    {
        printf("O sistema e indeterminado ou impossivel.");
        exit(0);
    }
    else
    {
        for(i=0; i<n; i++)
        {
            copiaMatriz(mat, aux, n);
            for(j=0; j<n; j++)
            {
                *(*(aux + j) + i) = *((double*)getItem(res, j));
            }
            d1 = determinante(aux, n);
            Dx = d1/D;
            printf("%c = %0.2lf \n", *((char*)getItem(vars, i)), Dx);
        }
    }
    free(res);
    free(vars);
}

/**
    Aloca na memória espaço para uma matriz nxn e a retorna
    @param n: a ordem da matriz
    @return mat: a matriz NxN já alocada
*/
double ** inicializaMatriz(int n)
{

    int i, j;
    double ** mat;
    mat = (double **)malloc(pow(n, 2)*sizeof(double));
    for(i=0; i<n; i++)
    {
        *(mat + i) = (double *) malloc(n*sizeof(double));
        for(j=0; j<n; j++)
        {
            *(*(mat+i) + j) = 0;
        }
    }

    return mat;
}

/**
    copia a matriz src para a matriz cpy por referência
    @param src: a matriz a ser copiada
    @param cpy: a cópia a ser gerada
*/
void copiaMatriz(double** src, double** cpy, int n)
{
    int i,j;
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            *(*(cpy + i) + j)= *(*(src+i) + j);
}

/**
    Calcula o determinante de uma matriz A, de ordem n
    @param A: a matriz cujo determinante deve ser calculado
    @param n: a ordem da matriz cujo determinante deve ser calculado
    @return det: o determinante da respectiva matriz
*/
double determinante(double **A, int n)
{
    double **Minor;
    int i,j,k,c1,c2;
    double det = 0;
    int O=1;

    Minor = inicializaMatriz(n);

    if(n==1)
    {
        return *(*(A));
    }
    if(n == 2)
    {
        det = *(*A)* (*(*(A+1) + 1))-(*(*A + 1))* *(*(A+1));
        return det;
    }
    else
    {
        for(i = 0 ; i < n ; i++)
        {
            c1 = 0, c2 = 0;
            for(j = 0 ; j < n ; j++)
            {
                for(k = 0 ; k < n ; k++)
                {
                    if(j != 0 && k != i)
                    {
                        *(*(Minor+c1) + c2)= *(*(A+j) + k);
                        c2++;
                        if(c2>n-2)
                        {
                            c1++;
                            c2=0;
                        }
                    }
                }
            }
            det = det + O*(*(*A + i)*determinante(Minor,n-1));
            O=-1*O;
        }
    }
    free(Minor);
    return det;
}

/**
    insere o valor correspondente ao endereço *inf na lista correspondente ao endereço *lis
    @param lis: o ponteiro para a lista lis em que o valor deve ser inserido
    @param inf: o ponteiro para o valor que deve ser inserido
    @param sz: o tamanho do tipo do dado a ser inserido na lista
    @return a posição em que foi inserido o valor
*/
char insira (Lista* lis, void* inf, size_t sz)
{
    char c;
    int i = 1;
    if (lis->inicio==NULL)
    {
        lis->inicio       = (No*)malloc(sizeof(No));
        lis->inicio->info = malloc(sz);
        memmove(lis->inicio->info, inf, sz);
        lis->inicio->prox = NULL;
        return i-1;
    }
    No* atual = lis->inicio;
    c = *((char*) lis->inicio->info);
    while(atual->prox != NULL)
    {
        i++;
        atual = atual->prox;
        c = *((char*) atual->info);
    }

    atual->prox = (No*)malloc(sizeof(No));
    atual->prox->info= malloc(sz);
    memmove(atual->prox->info, inf, sz);
    atual->prox->prox= NULL;
    free(inf);
    return i;
}

/**
    busca o valor armazenado no indice especificado a lista com endereço correspondente ao especificado
    @param lis: o ponteiro para a lista aonde deve-se buscar o valor
    @param index: o indice correspondente ao dado solicitado na lista
    @return o endereço correspondente ao valor armazenado na posição index da lista lis
*/
void* getItem(Lista *lis, int index)
{
    int i;
    No* point = lis->inicio;
    for(i=0; i<index; i++)
    {

        point = point->prox;
    }
    return point->info;
}
