# include <cmath>
#include <algorithm>
#include <vector>
#include <gmp.h>
#include <math.h>
#include <cassert>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <random>
# include <ctime>
# include <cstring>

using namespace std;

# define N 29835457
// # define N 29835457
# define NVARS 10
# define PXOVER 0.95

# define PMUTATION 0.5
# define POPSIZE 50
# define MAXGENS 100

struct genotype
{
  vector<int> gene;
  int fitness;
  double rfitness;
  double cfitness;

  genotype() {
    fitness = 864;
    rfitness = 0.0;
    cfitness = 0.0;
  }
};

vector<struct genotype> population;
vector<struct genotype> newpopulation;

vector<int> tobit(int n);
void initialize ();
int tobeten(vector<int> &x);
void evaluate();
void keep_the_best();
bool checkpropety(int x);
void Xover ( int one, int two, long &seed );
vector<int> made_crom_with_prop(vector<int> num);
void mutate ( long &seed );
void selector (long &seed );
void crossover ( long &seed );
void elitist ( );
void printpopulation();
bool sort_by_fitness( const genotype & lhs, const genotype & rhs );
bool compare_by_fitness( const genotype & lhs, const genotype & rhs );
void check();
bool checksolution();

// вид хромосомы - перевод из 10 СИ в 2 СИ
// возвращает вектор битов
vector<int> tobit(int n){
vector <int> vec;
    while (n)
    {
        vec.push_back(n%2);
        n /= 2;
    }
    reverse(vec.begin(), vec.end());
    return vec;
}

// изначальное заполнение популяции хромосомами
void initialize ()
{
  int i;
  ifstream input;
  int j;
  vector<int> num;
  vector<int> Nbit = tobit(N);
  cout << "Number fo divided"  << endl;
  int size_chrom = Nbit.size()/2;
  cout << "Lenght of chromosome = " << Nbit.size()/2 << endl;
  for (int j = 0 ; j < Nbit.size() ; j++){
      cout << Nbit[j];
  }
 cout << endl;
 cout << "first_chrom" << endl;
 vector<int> first_chrom (size_chrom, 0);
 first_chrom[0] = 1;
   for (int j = 0 ; j < first_chrom.size() ; j++){
      cout << first_chrom[j];
  }
   cout << endl;
   cout << "first chrome in ten" << tobeten(first_chrom) << endl;
  // i = tobeten(first_chrom)-1;
  for (i = 3 ; i < sqrt(N) ; i++){
    num.push_back(6 * i + 1);
    num.push_back(6 * i - 1);
  }
    for ( j = 2; j < POPSIZE+2; j++ )
    {
      genotype gt;

      // cerr << population.size() << endl;
      gt.fitness = 0;
      gt.rfitness = 0;
      gt.cfitness = 0;
      gt.gene = tobit(num[j]);
      population.push_back(gt);
    }
  cout << "population.size " << population.size() << endl;
  return;
}

int tobeten(vector<int> &x){
  int myNumInInteger = 0 ;
      for (int i = x.size(); i> 0; i--) {
        // cout << "x[i] = " << x[i];
        // cout << "pow(2, i)" << pow(2, i-1) << "i " << x.size()-i << "x[i] " << x[x.size()-i] << endl;
        myNumInInteger = myNumInInteger + pow(2, i-1 )*x[x.size()-i];
    }
    return myNumInInteger;
}

// расчёт фитнесс-функции для каждой хромосомы в популяции
void evaluate ( )
{
  int gen;
  vector<int> x;

  for (int member = 0; member < population.size(); member++ )
  {
    int myNumInInteger = 0;
    x = population[member].gene;
    myNumInInteger = tobeten(x);
    // cout << "myNumInInteger" << myNumInInteger << endl;
    if (myNumInInteger != 0 ){
    population[member].fitness = N% myNumInInteger;
    }
    else { population.erase(population.begin() + member); }
    // cout << "population[member].fitness " << population[member].fitness << endl;
  }
  return;
}

// определение лучшей особи
void keep_the_best(){
  cout << "keep_the_best" << endl;
  int cur_best;
  int mem;
  int i;

  cur_best = 0;

  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    // cout << "POPSIZE" << POPSIZE << endl;
    // cout << "population[POPSIZE].fitness" << population[POPSIZE].fitness << endl;
    if ( population[POPSIZE].fitness > population[mem].fitness )
    {
      cur_best = mem;
      population[POPSIZE].fitness = population[mem].fitness;
    }
  }
}

double r8_uniform_ab ( double a, double b, long &seed )
{
  int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    // cerr << "\n";
    // cerr << "R8_UNIFORM_AB - Fatal error!\n";
    // cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}

int i4_uniform_ab ( int a, int b, long &seed )
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    // cerr << "\n";
    // cerr << "I4_UNIFORM_AB - Fatal error!\n";
    // cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}

bool checkpropety(int x){
  // x = 6*m+1
  // m = x -1 / 6
  int m;
  m = x - 1;
  if (m % 6 ==0){
    return 1;
  }
  return 0;
}

void printvector(vector<int> a){
  for (int i = 0; i<a.size(); i++){
      cout << a[i];
  }
  cout << endl;
}

void Xover ( int one, int two, long &seed )
{
  int i, t;
  int point;
  double firstpart, twopart;
  // cout << "1 " << population[one].gene.size() << endl;
  // cout << tobeten(population[one].gene) << endl;
  // cout << "2 " << population[two].gene.size() << endl;
  // cout << tobeten(population[two].gene) << endl;
  // cout << "ffffffffffff" << endl;
  // printvector(population[two].gene);
  int diff = newpopulation[one].gene.size() - newpopulation[two].gene.size();
  // cout << "DIFF " << diff << endl;
  if (diff == 0 ){
    point = newpopulation[one].gene.size() /2;
  }
  else if (abs(diff)> 100){
    cout << "Что-то пошло не так...." << endl;

    // cout << "GEN " << tobeten(population[two].gene) << endl;
    // cout << tobeten(population[one].gene) << endl;
  }
  else if (diff > 0){
    for (int k; k <= diff; k++){

    newpopulation[two].gene.insert(newpopulation[two].gene.begin(), 0);
    }
  } 
  else {
    // cout << "ЗАШЕЛ" << endl;
    diff = (-1)*diff;
    // cout << "DIFF " << diff << endl;
    for (int k=0; k <= diff; k++){
      newpopulation[one].gene.insert(newpopulation[one].gene.begin(), 0);
    }
    // cout << "SIZE1" << population[one].gene.size() << endl;
    // cout << "SIZE2" << population[two].gene.size() << endl;
  }
  

  // cout << "point" << point << endl;
//
//  Swap genes in positions 0 through POINT-1.
//
genotype popmutateone;
genotype popmutatetwo;
popmutateone.gene = newpopulation[one].gene;
popmutatetwo.gene = newpopulation[two].gene;
// cout << " population[one].gene[i " <<  population[one].gene.size() << endl;
// cout << " population[two].gene[i " <<  population[two].gene.size() << endl;
  for ( i = 0; i < point; i++ )
  {
    t = popmutateone.gene[i];
    popmutateone.gene[i] = popmutatetwo.gene[i];
    popmutatetwo.gene[i] = t;
    vector<int> num1 = popmutateone.gene;
    popmutateone.gene = made_crom_with_prop(num1);
    vector<int> num2 = popmutatetwo.gene;
    popmutatetwo.gene = made_crom_with_prop(num2);
  }
  // if (checkpropety(tobeten(popmutateone.gene))){
  //   population.push_back(popmutatetwo);
  // }
  // else if (checkpropety(tobeten(popmutatetwo.gene))){
  //    population.push_back(popmutateone);
  // }
  population.push_back(popmutatetwo);
  // cout << "FIRST" << tobeten(popmutateone.gene) << endl;
  population.push_back(popmutateone);
  // cout << "SECOND" << tobeten(popmutatetwo.gene) << endl;
  return;
}

vector<int> made_crom_with_prop(vector<int> num){
  vector<int> a;
    if (checkpropety(tobeten(num)) ){
        return num;
     }
    else {
      if((tobeten(num)/ 6 )==0){
        int number = rand() % (N/6);
        return tobit( number* 6 + 1);
      }
      else{
        a.push_back((tobeten(num)/ 6 ) * 6 + 1);
        a.push_back((tobeten(num)/ 6 ) * 6 - 1);
        int random_number = rand() % 2;
        return tobit(a[random_number]);
      }
    }
}

void mutate ( long &seed )
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;

  for ( i = 0; i < population.size(); i++ )
  {
    for ( j = 0; j < population[i].gene.size(); j++ )
    {
      x = r8_uniform_ab ( a, b, seed );
      // cout << "mutate      x      " << x << endl;
      if ( x < PMUTATION )
      {
      // cout << "ЗАШЕЛ" << endl;
      int last_value = population[i].gene.size()-1;
      int random_number = 1 + rand() % last_value;
      // cout << "random_number " << random_number << endl;
      vector<int> copy  = population[i].gene;
      copy[random_number] ^= 1;
      
      population[i].gene = made_crom_with_prop(copy);
      // cout << "NEW AFTER MUTATE   " << tobeten(population[i].gene) << endl;
      }
    }
  }

  return;
}
// надо поправить параметр p
void selector ( long &seed )
{
  vector<genotype> newpop;
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  int sum;
//
//  Find the total fitness of the population.
//
  sum = 0;
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    population[mem].rfitness = population[mem].fitness / sum;
  }
// 
//  Calculate the cumulative fitness.
//
  population[0].cfitness = population[0].rfitness;
  for ( mem = 1; mem < POPSIZE; mem++ )
  {
    population[mem].cfitness = population[mem-1].cfitness +       
      population[mem].rfitness;
  }
// 
//  Select survivors using cumulative fitness. 
//
  for ( i = 0; i < population.size(); i++ )
  { 
    // p = 23;
    p = r8_uniform_ab(a,b,seed);
    // cout << "p " << p << endl;
    if ( p > 0.5 )
    {
      newpopulation.push_back(population[i]);    
    }
  }
// 
//  Overwrite the old population with the new one.
//
// cout << "newpop.size()" << newpopulation.size() << endl;
  // for ( i = 0; i < newpopulation.size(); i++ )
  // {

  //   // cout << "new_pop_fit" << newpopulation[i].fitness << endl; 
  // }

  // cout << "chromosome for mutate" << population.size() << endl;
  return;     
}

void crossover ( long &seed )
{
  const double a = 0.0;
  const double b = 1.0;
  int mem;
  int one;
  int first = 0;
  double x;

  for ( mem = 0; mem < newpopulation.size(); ++mem )
  {
    x = r8_uniform_ab ( a, b, seed );
    // cout << "x         crossover             " << x << endl;
    if ( x < PXOVER )
    {
      // cout << "ЗАШЕЛ" << endl;
      ++first;

      if ( first % 2 == 0 )
      {
        // cout << "one " << one << "mem " << mem << endl; 
        Xover ( one, mem, seed );
      }
      else
      {
        one = mem;
      }
    // printpopulation();
    }
  }
  return;
}

void elitist ( )
{
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;

  best = population[0].fitness;
  worst = population[0].fitness;

  for ( i = 0; i < POPSIZE - 1; ++i )
  {
    if ( population[i+1].fitness < population[i].fitness )
    {

      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }

    }

  }
  if ( population[POPSIZE].fitness <= best )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[POPSIZE].gene[i] = population[best_mem].gene[i];
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
  }
  else
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[worst_mem].gene[i] = population[POPSIZE].gene[i];
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  } 

  return;
}

void printpopulation(){
    for (int i = 0; i< population.size(); i++)
    {
      cout << i << " genes value " << " in ten " << tobeten(population[i].gene) << " in bit ";
      for (int j = 0 ; j < population[i].gene.size(); j++)
      {
      cout << population[i].gene[j];
      }
      
      cout << " fitness " << population[i].fitness << " cfi" << endl;
    }
}

bool sort_by_fitness( const genotype & lhs, const genotype & rhs )
{
   return lhs.fitness < rhs.fitness;
}

bool compare_by_fitness( const genotype & lhs, const genotype & rhs )
{
   return lhs.fitness == rhs.fitness;
}

void check(){
  for (int i = 0; i < population.size(); i++){
    if (tobeten(population[i].gene) == 0 || tobeten(population[i].gene) == 1){
      population.erase(population.begin() + i);
    }
  }
}

bool checksolution(){
   for (int i = 0; i < population.size(); i++){
   if (population[i].fitness == 0){
     return 1;
   }

   }
   return 0;
}

struct PRNG
{
    unsigned seed = 0;
};

void initGenerator(PRNG& generator)
{
    // Получаем случайное зерно последовательности
    generator.seed = unsigned(std::time(nullptr));
}

// Генерирует число на отрезке [minValue, maxValue].
unsigned random(PRNG& generator, unsigned minValue, unsigned maxValue)
{
    // Проверяем корректность аргументов
    assert(minValue < maxValue);
    // Итеративно изменяем текущее число в генераторе
    generator.seed = (generator.seed * 73129 + 95121);
    
    // Приводим число к отрезку [minValue, maxValue]
    return (generator.seed % (maxValue + 1 - minValue)) + minValue;
}

void add_new_chrome(){
  int size_pop = population.size();
  PRNG generator;
  initGenerator(generator);
  if (size_pop != POPSIZE){
    if (size_pop > POPSIZE){
      population.erase(population.begin() + POPSIZE, population.end());
    }
    else{
      while (population.size() != POPSIZE)
      {
      cout << "add_new_chrome" << endl;
      unsigned t = random(generator, POPSIZE , sqrt(N))  ;
      cout << t << endl;
      genotype gt;
      gt.fitness = N% t;
      gt.rfitness = 0;
      gt.cfitness = 0;
      gt.gene = tobit(t);
      population.push_back(gt); 
      t = 0;
      }        
    }
  }
}
int main(){
    string filename = "input.txt";
    int generation;
    int i;
    long seed = N;
    initialize();
    evaluate();
    for ( generation = 0; generation < MAXGENS; generation++ )
  {
    // timestamp( );
    seed = 123456789;
    // initialize ( filename, seed );
    evaluate();
    cout << "__initialization__" << endl;
    printpopulation();
    if (generation == 0){
      check();
      if (checksolution() ==1)
      {
        cout << "cancel!!!!!!" << endl;
        break;
      }
    }
    keep_the_best();
    cout << "__keep_the_best__" << endl;
    // printpopulation();
    selector ( seed );
    cout << "__selector__" << endl;
    // printpopulation();
    crossover ( seed );
    cout << "__crossover__" << endl;
    // printpopulation();
    mutate (seed );
    cout << "__mutate__" << endl;
    // printpopulation();
    evaluate ( );
    cout << "__sort__" << endl;
    sort(population.begin(), population.end(), sort_by_fitness);
    // printpopulation();
    cout << "__compare__" << endl;
    // printpopulation();
    population.erase(std::unique(population.begin(), population.end(), compare_by_fitness), population.end());
    add_new_chrome();
    // printpopulation();
    check();
    cout << "N = " << N << endl;
    if (checksolution() ==1)
    {
      printpopulation();
      int k = tobeten(population[0].gene);
      cout << "gen" << generation << " cancel!!!!!!" << k << endl;
      cout << "N = " << N << "; 1 = " << k << "; N%k = " << N%k << "; 2 = " << N/k << endl; 
      break;
    }
    cout << generation << "  SIZE  " << population.size() << endl;
  }
    return 0;
}


