typedef
    struct No
    {
        void* info;
        struct No* prox;
    }No;

typedef
    struct Lista
    {
        No* inicio;
    }Lista;
Lista* res;/*a �ltima coluna do sistema (a dos dados ap�s o '=') que dever� ser substitu�da na matriz na hora do c�lculo*/
Lista* vars;/*as vari�veis do sistema, �teis tanto para inserir os valores na matriz, tanto para escrever o valor das vari�veis no fim do c�lculo*/

double determinante(double **, int);
double** inicializaMatriz(int);
char insira (Lista*, void*, size_t);
void* getItem(Lista*, int);
void resolve(double**, int);
