#include<stdio.h>
#include<stdlib.h>
#include <sys/time.h>
#include<math.h>
#include<time.h>


// Structure de nombre complexe
typedef struct complex_ {
    double real;
    double imag;
}complex;

// Structure de polynome a coefficients complexes
typedef struct polynomeComplex_ {
    int degree;
    complex* tab;
} polynomeComplex;

// Retourne true si n>0 et est une puissance de deux sinon retourne false
int estPuissanceDeDeux(int n){
    return (n > 0) && ((n & (n -1)) == 0);
}

// Retourne la plus petite puissance de deux x tel que x > n
int puissanceDeDeuxSuperieur(int n){
    int puissance = 1;
    while(puissance < n) {
        puissance <<= 1;
    }
    return puissance;
}

/*Arithmetique des nombres complexe  */

// addition de deux nombres complexe
//(a+ib) + (c+id) = ((a+c)+i(b+d)) 
complex add(complex* a, complex* b){
    complex result;
    result.real = a->real + b->real;
    result.imag = a->imag + b->imag;
    return result;
}

// soustraction de deux nombres complexe
//(a+ib) - (c+id) = ((a-c)+i(b-d)) 
complex sub(complex* a, complex* b){
    complex result;
    result.real = a->real - b->real;
    result.imag = a->imag - b->imag;
    return result;
}

// conjuguee d'un nombre complexe
//conjugate(a+ib) = a-ib
complex conjugate(complex *a){
    complex result;
    result.real= a->real;
    result.imag = -1 * a->imag;
    return result;
}

// multiplication de deux nombres complexe 
//(a+ib)  * (c+id) = ac - bd + i(ad + bc)
complex c_mult(complex *a, complex *b){
    complex result;

    double real_part = a->real * b->real - a->imag * b->imag;

    double imag_part = a->real * b->imag + a->imag * b->real;

    result.real = real_part;
    result.imag = imag_part;

    return result;
}
// Division nombre complexe par un entier
//(a+ib)/n = ((a/n)+i(b/n))
complex c_div(complex *a, int n){
    complex result;
    result.real = a->real/n;
    result.imag = a->imag/n;
    return result;
}

// fonction retournant le temps instant t
double wtime(){
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

// Affiche le polynome complexe 
void affichagePolynomeComplex(polynomeComplex* polynome){
    for(int i = 0; i< polynome->degree+1; i++){
        if(i==0){
            printf("(%f+i%f)X^%d ",polynome->tab[i].real,polynome->tab[i].imag,i);
        }
        else if(i == polynome->degree){
            printf("(%f+i%f)X^%d\n",polynome->tab[i].real,polynome->tab[i].imag,i);
        }else{
            printf("(%f+i%f)X^%d  ",polynome->tab[i].real,polynome->tab[i].imag,i);
        }
    }
}

// FFT utilisant un polynome entier en entrée
void FFT(int size, int* polynome, complex* result, int stride){
    if(size == 1){
        result[0].real = polynome[0];
        result[0].imag =0;
        return;
    }
    FFT(size/2,polynome,result, 2*stride);
    FFT(size/2,polynome + stride, result + size/2, 2*stride);
    complex p;
    complex q;
    complex w;
    w.imag = 0;
    w.real = 1;
    for(int j = 0; j< size/2; j++){
        p = result[j];
        q = c_mult(&result[j + size/2], &w);
        result[j] = add(&p,&q);// p + q
        result[j + size/2] = sub(&p,&q);// p - q
        w.real = cos((2*(j+1)*M_PI)/size);
        w.imag = -sin((2*(j+1)*M_PI)/size);
    }
}

// FFT utilisant un polynome complex en entrée
void FFTComplex(int size, complex* polynome, complex* result, int stride){
    if(size == 1){
        result[0] = polynome[0];
        return;
    }
    FFTComplex(size/2,polynome,result, 2*stride);
    FFTComplex(size/2,polynome + stride, result + size/2, 2*stride);
    complex p;
    complex q;
    complex w;
    w.imag = 0;
    w.real = 1;
    for(int j = 0; j< size/2; j++){
        p = result[j];
        q = c_mult(&result[j + size/2], &w);
        result[j] = add(&p,&q);// p + q
        result[j + size/2] = sub(&p,&q);// p - q
        w.real = cos((2*(j+1)*M_PI)/size);
        w.imag = -sin((2*(j+1)*M_PI)/size);
    }
}

// FFT utilisant un polynome entier en entree utilisant w déja calcule
void FFTOpti(int size, int* polynome, complex* result, int stride, complex** w){
     if(size == 1){
        result[0].real = polynome[0];
        result[0].imag =0;
        return;
    }
    FFTOpti(size/2,polynome,result, 2*stride,w+1);
    FFTOpti(size/2,polynome + stride, result + size/2, 2*stride,w+1);
    complex p;
    complex q;
    for(int j = 0; j< size/2; j++){
        p = result[j];
        q = c_mult(&result[j + size/2], &w[0][j]);
        result[j] = add(&p,&q);// p + q
        result[j + size/2] = sub(&p,&q);// p - q
    }
}

// FFT utilisant un polynome complex en entree utilisant w déja calcule
void FFTComplexOpti(int size, complex* polynome, complex* result, int stride, complex** w){
     if(size == 1){
        result[0] = polynome[0];
        return;
    }
    FFTComplexOpti(size/2,polynome,result, 2*stride,w+1);
    FFTComplexOpti(size/2,polynome + stride, result + size/2, 2*stride,w+1);
    complex p;
    complex q;
    for(int j = 0; j< size/2; j++){
        p = result[j];
        q = c_mult(&result[j + size/2], &w[0][j]);
        result[j] = add(&p,&q);// p + q
        result[j + size/2] = sub(&p,&q);// p - q
    }
}

// IFFT IMPLEMENTATION
polynomeComplex* IFFTApplication(complex* polynome, int degre){
    for(int i=0; i<degre+1; i++){
        polynome[i] = conjugate(&polynome[i]);
    }
    polynomeComplex* result = malloc(sizeof(polynomeComplex));
    result->degree = degre;
    result->tab = calloc(degre+1,sizeof(complex));
    FFTComplex(degre+1, polynome, result->tab,1);
    for(int i=0; i<degre+1; i++){
        result->tab[i] = c_div(&result->tab[i], degre+1);
    }

    return result;
}

// IFFT IMPLEMENTATION Optimisation pre-calcul de w
polynomeComplex* IFFTApplicationOpti(complex* polynome, int degre, complex** w){
    for(int i=0; i<degre+1; i++){
        polynome[i] = conjugate(&polynome[i]);
    }
    polynomeComplex* result = malloc(sizeof(polynomeComplex));
    result->degree = degre;
    result->tab = calloc(degre+1,sizeof(complex));
    FFTComplexOpti(degre+1, polynome, result->tab,1, w);
    for(int i=0; i<degre+1; i++){
        result->tab[i] = c_div(&result->tab[i], degre+1);
    }
    return result;
}

// Fonction verifie le polynome et execute la FFT sur le polynome
void FFTApplication(int* polynome, int degre){
    polynomeComplex* result = malloc(sizeof(polynomeComplex));
    printf("Test malloc\n");
    int size = degre+1;
    printf("degre = %d\n",size-1);

    // Verifie si la taille du polynome est une puissance de deux sinon ajoute des 0 pour completer
    if(!estPuissanceDeDeux(size)){
        
        printf("Pas une puissance de deux\n");
        size = puissanceDeDeuxSuperieur(size);
        printf("newSize = %d\n", size);

        polynome =realloc(polynome,size * sizeof(int));
        for(int i = degre+1; i < size; i++){
            polynome[i]= 0;
        }
        degre = size-1;
    }

    result->degree = degre;
    result->tab = calloc(size,sizeof(complex));
    FFT(size, polynome, result->tab, 1);
    affichagePolynomeComplex(result);
    IFFTApplication(result->tab,result->degree);

}

// Affiche le polynome complex selon le degre max
void affichagePolynome(polynomeComplex* polynome, int degreMax){
    for(int i = degreMax; i >=0; i--){
        if(i== degreMax){
            printf("[%d, ",(int) round(polynome->tab[i].real));
        }
        else if(i == 0){
            printf("%d]\n",(int) round(polynome->tab[i].real));
        }else{
            printf("%d, ",(int) round(polynome->tab[i].real));
        }
    }
}

// Affiche le polynome entier selon le degre max
void affichagePolynomeInt(int* polynome, int degreMax){
    for(int i = degreMax; i >=0; i--){
        if(i== degreMax){
            printf("[%d, ",polynome[i]);
        }
        else if(i == 0){
            printf("%d]\n",polynome[i]);
        }else{
            printf("%d, ",polynome[i]);
        }
    }
}

// multiplication de polynome utilisant la FFT
void multiplication(int* firstPolynome, int firstDegre, int* secondPolynome, int secondDegre){
    double temps = wtime();
    int n = firstDegre + secondDegre +1;
    int oldDegree = n-1;

    // Verification taille puissance de deux
    if(!estPuissanceDeDeux(n)){
        n = puissanceDeDeuxSuperieur(n);
    }

    firstPolynome =realloc(firstPolynome,n * sizeof(int));
    for(int i = firstDegre+1; i < n; i++){
        firstPolynome[i]= 0;
    }
    firstDegre = n -1;

    secondPolynome =realloc(secondPolynome,n * sizeof(int));
    for(int i = secondDegre+1; i < n; i++){
        secondPolynome[i]= 0;
    }
    secondDegre = n -1;

    // PRE-CALCUL DE w
    int tailletmp = n;
    int nombreDeDivisions = 0;
    int division[32]; //On suppose que le nombre de division est inferieur à 32

    // Calcul du nombre de division par 2 de la taille
    // Tant que l'entier n'est pas égal à 0
    while (tailletmp != 0) {

        // On stocke la division dans le tableau
        division[nombreDeDivisions] = tailletmp;

        // Effectuer une division par 2
        tailletmp = tailletmp >> 1;

        // Incrémenter le compteur de divisions
        nombreDeDivisions++;
    }

    // creation du tableau w
    complex** w = malloc(sizeof(complex*)*nombreDeDivisions);
    
    // Calcul de w pour la taille se trouvant dans division[i]
    for(int i = 0; i < nombreDeDivisions; i++){
        int divisionTmp = division[i];
        w[i] = malloc(sizeof(complex)*(divisionTmp/2));
        w[i][0].real = 1;
        w[i][0].imag = 0;
        for(int j = 1; j < (divisionTmp/2); j++){
            w[i][j].real = cos((2*j*M_PI)/divisionTmp);
            w[i][j].imag = -sin((2*j*M_PI)/divisionTmp);
        }
    }
 
    // FFT
    // O(n*log(n))
    polynomeComplex* resultFirstPoly = malloc(sizeof(polynomeComplex));
    resultFirstPoly->degree = firstDegre;
    resultFirstPoly->tab = calloc(firstDegre+1,sizeof(complex));
    FFTOpti(firstDegre+1, firstPolynome, resultFirstPoly->tab, 1, w);

    polynomeComplex* resultSecondPoly = malloc(sizeof(polynomeComplex));
    resultSecondPoly->degree = secondDegre;
    resultSecondPoly->tab = calloc(secondDegre+1,sizeof(complex));
    FFTOpti(secondDegre+1, secondPolynome, resultSecondPoly->tab, 1, w);
    
    // MULTIPLICATION COEFF PAR COEFF
    // O(n)
    for(int i = 0; i< n;i++){
        resultFirstPoly->tab[i] = c_mult(&resultFirstPoly->tab[i],&resultSecondPoly->tab[i]);
    }

    // IFFT R
    // O(n*log(n))
    polynomeComplex* resultR;
    resultR = IFFTApplicationOpti(resultFirstPoly->tab,resultFirstPoly->degree, w);

    double endTime = wtime();

    affichagePolynome(resultR, oldDegree);
    printf("\nMultiplication FFT en %fs\n", endTime - temps);
}

// Multiplication naive sur des polynomes a coefficents entiers en O(n^2)
int* naiveMultiplication(int* firstVector, int firstDegre, int* secondVector, int secondDegre){
    double temps = wtime();
    int degree = firstDegre + secondDegre;
    int* result = malloc(sizeof(int) * degree +1);
    for(int i = 0;i <= firstDegre; i++){
        for(int j = 0;j <= secondDegre; j++){
            result[i+j] += firstVector[i] * secondVector[j];
        }
    }
    printf("Multiplication naive en %fs\n\n", wtime() - temps);
    return result;    
}
// genere un polynome entier 
void polynomegenerator(int* tab, int taille){
    for(int i=0; i<taille;i++){
        tab[i] = rand() % (100 - 0 +1) + 0;
    }

}
int main(int argc, char* argv[]){
    srand(time(NULL));

    // Choix de la taille du polynome
    // taille par défaut
    int taille = 11;

    if (argc >= 2) {
        if (atoi(argv[1]) >0){
            taille = atoi(argv[1]);
        }
    }
    
    printf("Generation des polynomes\n");

    int* first = malloc(sizeof(int)*taille);
    polynomegenerator(first, taille);


    int* second = malloc(sizeof(int)*taille);
    polynomegenerator(second, taille);

    printf("\nExecution de la méthode naive\n");
    int* naiveResult = naiveMultiplication(first, taille-1, second, taille-1);
    affichagePolynomeInt(naiveResult, (taille-1 + taille-1));

    printf("\nExecution de la méthode via la FFT\n\n");
    multiplication(first, taille-1, second, taille-1);

    return 0;
}