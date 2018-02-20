#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
//#define Imag csqrt(-1)
typedef struct{
    double complex **ptr;
    int row;
    int col;
}Matrix;

int getNullIndex(double complex *arr, int len);
void printMatrix(Matrix *m);
void emtyMatrix(Matrix *m);
void freeMatrix(Matrix *m);
Matrix* createMatrix(int row, int col);
Matrix* teilMatrix(Matrix*m, int row, int col);
void lambdaMatrix(Matrix*m, double complex lambda);
Matrix* copyMatrix(Matrix *m);
void swapRow(Matrix *m, int one, int two);
void sortByNullIndex(Matrix *m);
void subMultRow(Matrix *m, int one, double complex lone, int two, double complex ltwo);
void put(Matrix *m, int row, int col, double complex value);
void gauss(Matrix*m);
double complex det2x2(Matrix*m);

void error(char str[]){
    printf("[Error]%s\n",str);
    exit(-1);
}

int getNullIndex(double complex *arr, int len){
    double complex *p = arr;
    int cnt = 0;
    short pivot = 0;
    for (int i = 0; i < len; ++i, p++) {
        if(creal(*p) == 0.0 && cimag(*p) == 0.0)cnt++;
        else if(pivot == 0 ) pivot = 1;
        else return cnt;
    }
    return cnt;
}

void printMatrix(Matrix *m){
    printf("---\n");
    double complex **p = m->ptr;
    for (int i = 0; i < m->row; ++i) {
        double complex *t = *p++;
        printf("%p]\t",t);
        for (int j = 0; j < m->col; ++j) {
            double complex r = *t++;
            printf("%.0lf+%.0lf*i\t", creal(r), cimag(r));
        }
        printf("\n");
    }
}

void emtyMatrix(Matrix *m){
    for (int i = 0; i < m->row; ++i) {
        for (int j = 0; j < m->col; ++j){
            *(*(m->ptr+i)+j) = (double complex) 0.0;
        }
    }
}

void freeMatrix(Matrix *m){
    if(m == NULL)return;
    double complex **p = m->ptr;
    for (int i = 0; i < m->row; ++i, p++) {
        free(*p);
    }
    free(m->ptr);
    free(m);
}

Matrix* createMatrix(int row, int col){
    double complex **p = NULL;
    Matrix *m = malloc(sizeof(Matrix));
    if(m == NULL)error("createMatrix");
    m->ptr = malloc(sizeof(double complex *)*row);
    p = m->ptr;
    if(m->ptr == NULL)exit(-1);
    for (int i = 0; i < row; ++i, p++) {
        *p = malloc(sizeof(double complex)*col);
        if(m->ptr+i == NULL)exit(-1);
    }
    m->col = col;
    m->row = row;
    emtyMatrix(m);
    return m;
}

Matrix* teilMatrix(Matrix*m, int row, int col){
    Matrix *n = createMatrix(m->row-1, m->col-1);
    for (int i = 0; i < m->row; ++i) {
        if(i == row)continue;
        for (int j = 0; j < m->col; ++j) {
            if(j == col)continue;
            if(i<row && j<col)
                *(*(n->ptr+i)+j) = *(*(m->ptr+i)+j);
            else if (i<row && j>col)
                *(*(n->ptr+i)+j-1) = *(*(m->ptr+i)+j);
            else if(i>row && j < col)
                *(*(n->ptr+i-1)+j) = *(*(m->ptr+i)+j);
            else
                *(*(n->ptr+i-1)+j-1) = *(*(m->ptr+i)+j);
        }
    }
    return n;
}

double complex laplace(Matrix*m){
    Matrix *x = NULL;
    double complex det = 0.0;
    if(m->row!=m->col)error("laplace");
    if(m->col<=2)return det2x2(m);
    for (int i = 0; i < m->row; ++i) {
        x = teilMatrix(m,i,0);
        det += **(m->ptr+i)* cpow(-1,i)* laplace(x);
        printf("<%.0lf%.0lf\n", creal(det),cimag(det));
        freeMatrix(x);
    }
    return det;
}


double complex det2x2(Matrix*m){
    double complex a,b,c,d,ad,bc;
    a = *(*(m->ptr));
    b = *(*(m->ptr)+1);
    c = *(*(m->ptr+1));
    d = *(*(m->ptr+1)+1);
    ad = a*d;
    bc = b*c;
  /*  printf("[%.0lf:%.0lf],[%.0lf:%.0lf],[%.0lf:%.0lf],[%.0lf:%.0lf]\n",
           creal(a),cimag(a), creal(b),cimag(b), creal(c), cimag(c), creal(d),cimag(d));
    printf(">%.0lfi%.0lf\n", creal(a*d),cimag(a*d));
    printf("<%.0lfi%.0lf\n", creal(b*c),cimag(b*c));
    printf("%.0lf,%.0lf\n", creal(ad-bc),cimag(ad-bc));
*/
     return ad-bc;
    //return *(*(m->ptr))**(*(m->ptr+1))-*(*(m->ptr+1)+1)**(*(m->ptr)+1)
}

void lambdaMatrix(Matrix*m, double complex lambda){
    for (int i = 0; i < m->col; ++i) {
        for (int j = 0; j < m->col; ++j) {
            *(*(m->ptr+i)+j) *= lambda;
        }
    }
}

Matrix* copyMatrix(Matrix *m){
    Matrix *t = createMatrix(m->row, m->col);
    for (int i = 0; i < m->row; ++i) {
        for (int j = 0; j < m->col; ++j) {
            *((*(t->ptr+i))+j) = *((*(m->ptr+i))+j);
        }
    }
    return t;
}

Matrix* mult(Matrix*one, Matrix*two){
    double complex **p,*a;
    Matrix*m=createMatrix(one->row, two->col);
    p = m->ptr;
    for (int i = 0; i < m->row; ++i) {
        for (int j = 0; j < m->col; ++j) {
            for (int k = 0; k < one->row; ++k) {
                *((*(p+i))+j) += (*((*(one->ptr+i))+k)) * *((*(two->ptr+k))+j);
            }
        }
    }
    return m;
}

void swapRow(Matrix *m, int one, int two){
    double complex *a = *(m->ptr+one);
    *(m->ptr+one) = *(m->ptr+two);
    *(m->ptr+two) = a;
}

// sort mit billigem sortierer TODO
void sortByNullIndex(Matrix *m){
    double complex **p = m->ptr;
    for (int i = 1; i < m->row; ++i) {
        if (getNullIndex(*(p+i), m->col) < getNullIndex(*(p+i-1), m->col)){
            swapRow(m, i, i-1);
            i=1;
        }
    }
}

void subMultRow(Matrix *m, int one, double complex lone, int two, double complex ltwo){
    double complex *a,*b;
    a = *(m->ptr+one);
    b = *(m->ptr+two);
    for (int i = 0; i < m->col; ++i) {
        *(a+i) = lone**(a+i) - ltwo* *(b+i);
    }
}

void put(Matrix *m, int row, int col, double complex value){
    *((*(m->ptr+row))+col) = value;
}

Matrix* transponiert(Matrix *m){
    Matrix *t = createMatrix(m->col, m->row);
    for (int i = 0; i < t->row; ++i) {
        for (int j = 0; j < t->col; ++j) {
            *(*(t->ptr + i) + j) = *(*(m->ptr + j) + i);
        }
    }
    return t;
}

void gauss(Matrix*m){
    //Matrix *b = copyMatrix(m);
    //emtyMatrix(b);
    double complex a,b;
    for (int c = 0; c < m->col; ++c) {
        for (int r = c+1; r < m->row; ++r) {
            //double complex temp = *(*(m->ptr+)+);
            //addMultRow(m, r, *(*(m->ptr+r-1)+c), r-1, *(*(m->ptr+r)+c));
            a = *(*(m->ptr+r-1)+c);
            b = *(*(m->ptr+r)+c);
            if(a != 0.0 && b != 0.0){
                subMultRow(m, r, a, r-1, b);
            }
            //printMatrix(m);
        }
    }
    sortByNullIndex(m);
}

int main(int argc, char* argv[]) {
    Matrix *t = createMatrix(3,3);
    Matrix *q = createMatrix(2,2);
    put(q,0,0,-I);
    put(q,0,1,2);
    put(q,1,0,3);
    put(q,1,1,4);

    printMatrix(q);

    double complex d = det2x2(q);
    printf("det=%.0lf,%.0lf\n",creal(d), cimag(d));

    put(t, 0, 0, 5.0 );
    put(t, 0, 1, 6.0 );
    put(t, 0, 2, 7.0 );

    put(t, 1, 0, 1.0);
    put(t, 1, 1, 2.0);
    put(t, 1, 2, 3.0);

    put(t, 2, 2, 1.0);

    Matrix * m = transponiert(t);
    //printMatrix(t);

    //printf("det=%lf\n", laplace(t));
    //Matrix *x = teilMatrix(t,1,1);
    //printMatrix(x);
    //printf("---\n");
    //zeilenStufenForm(t);
    //gauss(t);
    //printMatrix(t);
    //rowAddMultRow(t, 0, 0.5, 0);
    //multRowByRow(t, 2, 3.5, 0 );
    //printMatrix(t);
    //m = copyMatrix(t);
    //multRowByLambda(t, 0, 2.0);
    //printMatrix(t);
    //t = mult(t,t);
    //printMatrix(m);
    //printMatrix(t);
    freeMatrix(t);
    //freeMatrix(x);
    printf("finish");
    return 0;
}