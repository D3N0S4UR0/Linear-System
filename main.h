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
Lista* res;/*a última coluna do sistema (a dos dados após o '=') que deverá ser substituída na matriz na hora do cálculo*/
Lista* vars;/*as variáveis do sistema, úteis tanto para inserir os valores na matriz, tanto para escrever o valor das variáveis no fim do cálculo*/

double determinante(double **, int);
double** inicializaMatriz(int);
char insira (Lista*, void*, size_t);
void* getItem(Lista*, int);
void resolve(double**, int);
