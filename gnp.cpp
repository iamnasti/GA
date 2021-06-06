#include <stdio.h>
# include <iostream>
#include <gmp.h>
#include <vector>
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

#include <memory>

using namespace std;

# define POPSIZE 2000
# define MAXGENS 2000
# define PCROSSOVER 1
# define PXOVER 1
# define PMUTATION 1

class BigInt {
public:
  BigInt() {
    mpz_init(v);
    int i = 0;
  }

  BigInt(const BigInt &num) {
    // mpz_init(v);
    mpz_init_set(v, num.v);
  }

  BigInt(const string &s) {
    mpz_init(v);
    mpz_t base, num;
    mpz_init_set_ui(base, 10);
    mpz_init(num);

    for (char c : s) {
      mpz_mul(v, v, base);
      mpz_set_ui(num, c-'0');
      mpz_add(v, v, num);
    }

    mpz_clear(base);
    mpz_clear(num);
  }

  ~BigInt(){
    mpz_clear(v);
  }

  // BigInt operator+(const mpz_t& b) const {
  //   BigInt res;
  //   mpz_add(res.v, v, b);
  //   return res;
  // }

  BigInt operator+(const BigInt& b) const {
    BigInt res;
    mpz_add(res.v, v, b.v);
    return res;
  }

  BigInt operator+(const int b) const {
    BigInt res;
    BigInt num;
    num = b;
    mpz_add(res.v, v, num.v);
    return res;
  }

  BigInt operator-(const BigInt& b) const {
    BigInt res;
    mpz_sub(res.v, v, b.v);
    return res;
  }

  BigInt pow(const unsigned int& b) const {
    BigInt res;
    mpz_pow_ui(res.v, v, b);
    return res;
  }

  BigInt root(const unsigned int& b) const {
    BigInt res;
    mpz_root(res.v,v, b);
    return res;
  }

  BigInt xor_num(const BigInt& b) const {
    BigInt res;
    mpz_xor(res.v, v, b.v);
    return res;
  }

  auto is_prime() const {
    return mpz_probab_prime_p(v, 40);
  }

  BigInt random_num() const{
    BigInt res;
    mpz_random (res.v, 500);
    return res;
  }

  BigInt operator-(const int b) const {
    BigInt res;
    BigInt num;
    num = b;
    mpz_sub(res.v, v, num.v);
    return res;
  }

  BigInt operator*(const BigInt& b) const {
    BigInt res;
    mpz_mul(res.v, v, b.v);
    return res;
  }

  BigInt operator*(const int b) const {
    BigInt res;
    BigInt num;
    num = b;
    mpz_mul(res.v, v, num.v);
    return res;
  }

  const BigInt operator/(const BigInt& b) const {
    BigInt res;
    mpz_fdiv_q(res.v, v, b.v);
    return res;
  }

  const BigInt operator/(const int b) const {
    BigInt res;
    BigInt num;
    num = b;
    mpz_fdiv_q(res.v, v, num.v);
    return res;
  }

  const BigInt operator%(const int b) const {
    BigInt res;
    BigInt num;
    num = b;
    mpz_mod(res.v, v, num.v);
    return res;
  }

  const BigInt operator%(const BigInt& b) const {
    BigInt res;
    mpz_mod(res.v, v, b.v);
    return res;
  }

  const BigInt operator/=(const BigInt& b) {
    BigInt res;
    mpz_fdiv_q(v, v, b.v);
    return *this;
  }

  const BigInt operator/=(const unsigned int b) {
    BigInt num;
    num = b;
    BigInt res;
    mpz_fdiv_q(v, v, num.v);
    return *this;
  }

  BigInt& operator=(const BigInt& b) {
      mpz_set(v, b.v);
      return *this;
  }

  BigInt& operator=(const unsigned int b) {
      mpz_set_ui(v, b);
      return *this;
  }

  // BigInt& operator=(const mpz_t& b) {
  //     mpz_set(v,b);
  //     return *this;
  // }
  bool operator==(const BigInt& b) const {      
      return mpz_cmp(v,b.v) == 0;
  }

  bool operator>(const BigInt& b) const  {      
    if (mpz_cmp(v,b.v) == 1){
      return true;
    }
    else {return false;}
  }

  bool operator<(const BigInt& b) const {      
    if (mpz_cmp(v,b.v) == -1){
      return true;
    }
    else {return false;}
  }




  friend std::ostream& operator<<(std::ostream &os, const BigInt &b) {
    char s[256];
    mpz_get_str(s,10,b.v);
    os << s;
    return os;
  }


private:
  mpz_t v;

};

struct genotype
{
  vector<BigInt> gene;
  BigInt fitness;
  BigInt rfitness;
  BigInt cfitness;

  genotype() {
    fitness = 864;
    rfitness = 0.0;
    cfitness = 0.0;
  }
};

vector<struct genotype> population;
vector<struct genotype> newpopulation;
vector <BigInt> tobit(BigInt n);
void initialize (BigInt N);
BigInt tobeten(vector<BigInt> &x);
void evaluate (BigInt N );
void printvector(vector<BigInt> a);
void check();
int checksolution();
bool sort_by_fitness( const genotype & lhs, const genotype & rhs );
void selector ();
void crossover (BigInt N);
void Xover ( int one, int two, BigInt N);
vector<BigInt> made_crom_with_prop(vector<BigInt> num, BigInt N);
bool checkpropety(BigInt x);
bool compare_by_fitness( const genotype & lhs, const genotype & rhs );
void add_new_chrome(BigInt N);

void GeneSizes() {
  for (int i = 0; i < population.size(); i++) {
  cout << population[i].gene.size() << " ";
    if (population[i].gene.size() > 1000) {
      cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ " ;
    }
  }
  cout << "КОНЕЦ" << endl;
  cout << endl;
}

void NewGeneSizes() {
  for (int i = 0; i < newpopulation.size(); i++) {
  cout << newpopulation[i].gene.size() << " ";
    if (newpopulation[i].gene.size() > 1000) {
      cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ " ;
    }
  }
  cout << endl;
}

vector <BigInt> tobit(BigInt n){
  // gmp_printf("N : %Zd", n.getValue());
  vector <BigInt> vec;
  BigInt zero;
  while (!(n == zero))
  {
      vec.push_back(n%2);
      // gmp_printf("_______%Zd : %d\n", n, vec.size());
      n /= 2;
  }
   reverse(vec.begin(), vec.end());
   if (vec.size() > 1000) {
    for (int i = 0; i < 200; i++) {
    }
  }
  
  return vec;
}

void initialize (BigInt N){
  int i;
  int j;
  vector<BigInt> num;
  vector<BigInt> Nbit = tobit(N);
  int size_chrom;
  BigInt zero;
  size_chrom = Nbit.size()/2;
  // cout << "Lenght of chromosome = " << Nbit.size()/2 << endl;
  // for (int j = 0 ; j < Nbit.size() ; j++){
  //     gmp_printf("%Zd", Nbit[j].getValue());
  // }
//  cout << endl;
//  cout << "first_chrom" << endl;
 vector<BigInt> first_chrom (size_chrom, zero);
 first_chrom[0] = 1;
 for (j = 0 ; j < first_chrom.size() ; j++){
      gmp_printf("%Zd", first_chrom[j]);
  }
   cout << endl;
  //  cout << "first chrome in ten" << tobeten(first_chrom) << endl;
 // i = tobeten(first_chrom)-1;
  BigInt i_1, max;
  unsigned int two;
  two = 2;
  // mpz_root(max_mp,N.getValue(), two);
  // gmp_printf("%Zd", max_mp);
  max = N.root(two);
  ifstream in("initial_number.txt");
  string line{};
  string number{};
   while (getline(in, line))
    {
        stringstream strStream(line);
        while (getline(strStream, number, ','))
        {
          BigInt ex;
          ex = atoi(number.c_str());
            num.push_back(ex);
        }

    }

  for ( j = 2; j < POPSIZE+2; j++ )
  {
    genotype gt;

//       // cerr << population.size() << endl;
      gt.fitness = 0;
      gt.rfitness = 0;
      gt.cfitness = 0;
      gt.gene = tobit(num[j]);
      population.push_back(gt);
    }
//   cout << "population.size " << population.size() << endl;
  return;
}

BigInt tobeten(vector<BigInt> &x){
  BigInt two, myNumInInteger, num1, pr,  num, res;
  two = 2;
  for (int i = x.size(); i> 0; i--) {
        // unsigned long int exp = i - 1;
        num1 = 0;
        pr = 0;
        res = two.pow(i-1);
        // mpz_pow_ui(res, two.getValue(), i-1 );
        num1 = res;
        pr = num1*x[x.size()-i];
        // cout << "x[i] = " << x[i];
        // cout << "pow(2, i)" << pow(2, i-1) << "i " << x.size()-i << "x[i] " << x[x.size()-i] << endl;
        num = num + pr;
        // mpz_add(num,num, pr.getValue());
    }
  myNumInInteger = num;
    return myNumInInteger;
}


// расчёт фитнесс-функции для каждой хромосомы в популяции
void evaluate (BigInt N )
{
  int gen;
  vector<BigInt> x;

  for (int member = 0; member < population.size(); member++ )
  {
    BigInt myNumInInteger;
    BigInt zero;
    x = population[member].gene;
    myNumInInteger = tobeten(x);
    // cout << "myNumInInteger" << myNumInInteger << endl;
    if (!(myNumInInteger == zero) ){
    population[member].fitness = N% myNumInInteger;
    }
    else { population.erase(population.begin() + member); }
    // cout << "population[member].fitness " << population[member].fitness << endl;
  }
  return;
}


void printvector(vector<BigInt> a){
  cout << "print vector" << endl;
  for (int i = 0; i<a.size(); i++){
      cout << a[i] << ", " << endl;
  }
  cout << endl;
}

void printpopulation() {
    cout << population.size() << endl;
    for (int i = 0; i< population.size(); i++)
    {

      cout << i << " genes value " << endl;
      // if (i == 13) {
      //   auto a = tobeten(population[i].gene).getValue();
      //   cout << endl;
      // }
      cout << "in ten " << tobeten(population[i].gene) << " in bit " << endl;
      // for (int j = 0 ; j < population[i].gene.size(); j++)
      // {
      //   cout << "i:" << i << " j:" << j << endl;
      //   if (i == 12 && j == 4) {
      //     cout << endl;
      //     auto w = population[i].gene[j];
      //     cout << w << endl;
      //   }
      // gmp_printf("%Zd",population[i].gene[j].getValue());
      // }
      
      cout << " fitness " << population[i].fitness << " cfi" << endl;
    }
}

void check(){
  BigInt zero, one;
  one = 1;
  for (int i = 0; i < population.size(); i++){

    if (tobeten(population[i].gene) == zero || tobeten(population[i].gene) == one){
      population.erase(population.begin() + i);
    }
  }
}

int checksolution(){
  BigInt zero;
  // cout << "____________check___________" << endl;
  for (int i = 0; i < population.size(); i++){
  if (population[i].fitness == zero){
    cout << "ВЫВОД: " << i << " : " << population[i].fitness << endl;
    return i;
  }
  }
  return -100;
}

void keep_the_best(){
  cout << "keep_the_best" << endl;
  int cur_best;
  int mem;
  int i;
  // printpopulation();
  cout << "___________________________________" << endl;
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
  // printpopulation();
}

bool sort_by_fitness( const genotype & lhs, const genotype & rhs )
{
   return lhs.fitness < rhs.fitness;
}

void selector ()
{
  vector<genotype> newpop;
  // srand(time(NULL));
  int i;
  int j;
  int mem;
  double p;
  BigInt sum;
  // for ( mem = 0; mem < POPSIZE; mem++ )
  // {
  //   sum = sum + population[mem].fitness;
  // }
  // for ( mem = 0; mem < POPSIZE; mem++ )
  // {
  //   population[mem].rfitness = population[mem].fitness / sum;
  // }
  // population[0].cfitness = population[0].rfitness;
  // for ( mem = 1; mem < POPSIZE; mem++ )
  // {
  //   population[mem].cfitness = population[mem-1].cfitness +       
  //     population[mem].rfitness;
  // }
  for ( i = 0; i < population.size(); i++ )
  { 
    p = rand()%100/100.0;
    // cout << "p " << p << endl;
    if ( p < PCROSSOVER)
    {
      newpopulation.push_back(population[i]);    
    }
  }
// cout << "number chromosome for mutate" << newpopulation.size() << endl;
  return;     
}

void crossover (BigInt N)
{
  // srand(time(NULL));
  int mem;
  int one;
  int first = 0;
  double x;
  // cout << "Размер популяции " <<   newpopulation.size() << endl;
  for ( mem = 0; mem < newpopulation.size(); ++mem )
  {
    x = rand()%100/100.0;
    // cout << "x         crossover             " << x << endl;
    if ( x < PXOVER )
    {
      // cout << "ЗАШЕЛ" << endl;
      ++first;

      if ( first % 2 == 0 )
      {
        // cout << " one " << one << "mem " << mem << endl; 
        // cout << "Xover начался" << endl;
        Xover ( one, mem, N);
        // cout << "Xover закончился" << endl;
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

void Xover ( int one, int two, BigInt N)
{
  // NewGeneSizes();
  int i;
  BigInt t;
  int point;
  BigInt zero;
  double firstpart, twopart;
  cout << "OLD" << endl;
  cout << "one" << tobeten(population[one].gene) << endl;
  cout << "two" << tobeten(population[two].gene) << endl;
  // cout << "1 " << population[one].gene.size() << endl;
  // cout << tobeten(population[one].gene) << endl;
  // cout << "2 " << population[two].gene.size() << endl;
  // cout << tobeten(population[two].gene) << endl;
  // cout << "ffffffffffff" << endl;
  // printvector(population[two].gene);
  int diff = newpopulation[one].gene.size() - newpopulation[two].gene.size();
  // cout << "DIFF " << diff << endl;
  // cout << "SIZE1 " << newpopulation[one].gene.size() << endl;
  // cout << "SIZE2 " << newpopulation[two].gene.size() << endl;
  if (diff == 0 ){
  }
  else if (abs(diff)> 100){
    cout << "Что-то пошло не так...." << endl;

    // cout << "GEN " << tobeten(population[two].gene) << endl;
    // cout << tobeten(population[one].gene) << endl;
  }
  else if (diff > 0){
    for (int k = 0; k <diff; k++){
// (gdb) p k <---------------------------------------------------------------
// $2 = -13216

    newpopulation[two].gene.insert(newpopulation[two].gene.begin(), zero);
    }
  } 
  else {
    // cout << "ЗАШЕЛ" << endl;
    diff = (-1)*diff;
    // cout << "DIFF " << diff << endl;
    for (int k=0; k <diff; k++){
      newpopulation[one].gene.insert(newpopulation[one].gene.begin(), zero);
    }
    // cout << "SIZE1" << population[one].gene.size() << endl;
    // cout << "SIZE2" << population[two].gene.size() << endl;
  }
  point = newpopulation[one].gene.size() /2;
  // cout << "SIZE1 " << newpopulation[one].gene.size() << endl;
  // cout << "SIZE2 " << newpopulation[two].gene.size() << endl;
  // cout << "point" << point << endl;
//
//  Swap genes in positions 0 through POINT-1.
//

genotype popmutateone;
genotype popmutatetwo;
popmutateone.gene = newpopulation[one].gene;
popmutatetwo.gene = newpopulation[two].gene;
// newpopulation.clear();
// NewGeneSizes();
// cout << "==========";
// GeneSizes();
// cout << "==========";
// cout << "point: " << point << endl;
// cout << " population[one].gene[i " <<  population[one].gene.size() << endl;
// cout << " population[two].gene[i " <<  population[two].gene.size() << endl;
for ( i = 0; i < point; i++ )
  {
    // cout << "номер точки: " << i << endl;
    t = popmutateone.gene[i];
    // cout << "last" << endl;
    // cout << "SIZE ONE: " << popmutateone.gene.size() << " " << popmutateone.gene[i]<< endl;
    // cout << "SIZE two: " << popmutatetwo.gene.size() << " " << popmutatetwo.gene[i] << endl;
    popmutateone.gene[i] = popmutatetwo.gene[i];
    // cout << "last" << endl;
    popmutatetwo.gene[i] = t;
    // cout << t << endl;
  }
  vector<BigInt> num1 = popmutateone.gene;
  popmutateone.gene = made_crom_with_prop(num1, N);
  vector<BigInt> num2 = popmutatetwo.gene;
  auto g =  made_crom_with_prop(num2, N);
  popmutatetwo.gene = g;
  // cout << "обмен генами закончен" << endl;
  // if (checkpropety(tobeten(popmutateone.gene))){
  //   population.push_back(popmutatetwo);
  // }
  // else if (checkpropety(tobeten(popmutatetwo.gene))){
  //    population.push_back(popmutateone);
  // }
  population.push_back(popmutatetwo);
  cout << "popmutateone : " << tobeten(popmutateone.gene) << endl;
  cout << "popmutatetwo : " << tobeten(popmutatetwo.gene) << endl; 
  // << "popmutatetwo : " <<  popmutatetwo.gene.size() << endl;
  // cout << "FIRST" << tobeten(popmutateone.gene) << endl;
  population.push_back(popmutateone);
  // cout << "SECOND" << tobeten(popmutatetwo.gene) << endl;
  return;
}

vector<BigInt> made_crom_with_prop(vector<BigInt> num, BigInt N){
  vector<BigInt> a;
  BigInt zero;
  // auto t = time(NULL);
  // cout << "time(NULL) : " << t;
  // srand(t);
    if (checkpropety(tobeten(num)) ){
        return num;
     }
    else {
      if((tobeten(num)/ 6 )==zero){
        BigInt randin;
        randin = rand();
        BigInt number;
        number= randin % (N/6);
        return tobit(number* 6 + 1);
      }
      else{
        a.push_back((tobeten(num)/ 6 ) * 6 + 1);
        a.push_back((tobeten(num)/ 6 ) * 6 - 1);
        int random_number;
        random_number = rand() % 2;
        return tobit(a[random_number]);
      }
    }
}

bool checkpropety(BigInt x){
  // x = 6*m+1
  // m = x -1 / 6
  BigInt m, zero;
  m = x - 1;
  if (m % 6 == zero && (m.is_prime() == 2 || m.is_prime() == 1 )){
    return 1;
  }
  return 0;
}

void mutate (BigInt N)
{
  int i;
  int j;
  int k = 0;
  BigInt one;
  one = 1;
  double x;
  // cout << "Рамер популяции " << population.size() << endl;
  for ( i = 0; i < population.size(); i++ )
  {
  x = rand()%100/100.0;
      // cout << "mutate      x      " << x << endl;
  if ( x < PMUTATION )
  {
    k++;
      // cout << "ЗАШЕЛ : " << k << " : " << population[i].gene.size() << " : " << population.size() << endl;
      int last_value = population[i].gene.size()-1;
      int random_number = 1 + rand() % last_value;
      // cout << "random_number " << random_number << endl;
      vector<BigInt> copy  = population[i].gene;
      // mpz_xor(res, copy[random_number].getValue(), one.getValue());
      copy[random_number] = copy[random_number].xor_num(one);
      population[i].gene = made_crom_with_prop(copy, N);
      // cout << "NEW AFTER MUTATE   " << tobeten(population[i].gene) << endl;
  }
  }
  // cout << "количество мутирующих " << k << endl;
  return;
}

bool compare_by_fitness( const genotype & lhs, const genotype & rhs )
{
   return lhs.fitness == rhs.fitness;
}

void add_new_chrome(BigInt N){
  gmp_randstate_t s;
  int size_pop = population.size();
  if (size_pop != POPSIZE){
    if (size_pop > POPSIZE){
      population.erase(population.begin() + POPSIZE, population.end());
    }
    else{
      while (population.size() != POPSIZE)
      {
      // cout << "add_new_chrome" << endl;
      BigInt t, zero, maxim, one;

      unsigned int two;
      two = 2;
      one = 1;
      // mpz_root(max_size,N.getValue(), two); 
      maxim = N.root(two);
      while (t == zero | t == one){

        // mpz_random (rop, 500);
        // t = rop;
        t = t.random_num()% maxim;
      }
      // cout << "t " << t << endl;
      genotype gt;
      gt.fitness = N% t;
      gt.rfitness = zero;
      gt.cfitness = zero;
      gt.gene = tobit(t);
      population.push_back(gt); 
      t = zero;
      }        
    }
  }
}


int main() {
  // auto t = time(NULL); // 1620675012
  // cout << "time(NULL) : " << t;
  srand(time(NULL));

  BigInt N("42336478013");
  
  int generation;
  // N = 6984871;
  // N = 10909343;
  // N = 29835457;
  // N = 392913607;
  // N = "5325280633";
  initialize(N);
  evaluate(N);
  cout << "__initialization__" << endl;
  // printpopulation();
  for ( generation = 0; generation <= MAXGENS; generation++ )
  {
    cout << "GENERATION: " << generation << endl;
    // cout << "POPILATION SIZE: " << population.size() << endl;
    // cout << "NEWPOPILATION SIZE: " << newpopulation.size() << endl;
    if (generation == 0){
    check();
    // printpopulation();
    int num_solution = checksolution();
      if (num_solution >= 0)
      {
        BigInt k;
      printpopulation();
      k = tobeten(population[num_solution].gene);
      cout << "result: gen" << generation << " cancel!!!!!!" << k << endl;
      cout << "result: N = " << N << "; 1 = " << k << "; N%k = " << N%k << "; 2 = " << N/k << endl; 
      break;
      }
    }
  
  // keep_the_best();
  // return 0;
  // cout << "__keep_the_best__: " << generation << endl;
  sort(population.begin(), population.end(), sort_by_fitness);
  // cout << "__keep_the_best__" << generation << endl;
  // printpopulation();
  // cout << "__selector__start" << endl;
  selector ();
  // cout << "__selector__" << endl;
  // printpopulation();
  crossover (N);
  // cout << "__crossover__" << endl;
  evaluate (N);
  sort(population.begin(), population.end(), sort_by_fitness);
  population.erase(std::unique(population.begin(), population.end(), compare_by_fitness), population.end());
  check();
  int num_solution = checksolution();
  if (num_solution >= 0)
  {
    BigInt k;
    printpopulation();
    k = tobeten(population[num_solution].gene);
    cout << "result: gen" << generation << " cancel!!!!!!" << k << endl;
    cout << "result: N = " << N << "; 1 = " << k << "; N%k = " << N%k << "; 2 = " << N/k << endl; 
    break;
  }
  mutate (N);
  // printpopulation();
  // cout << "__mutate__" << endl;
  evaluate (N);
  // cout << "__sort__" << endl;
  sort(population.begin(), population.end(), sort_by_fitness);
  population.erase(std::unique(population.begin(), population.end(), compare_by_fitness), population.end());
  // printpopulation();
  add_new_chrome(N);
  // printpopulation();
  check();
  cout << "NEWPOPILATION SIZE: " << newpopulation.size() << endl;
  newpopulation.clear();
  cout << "N = " << N << endl;
  int num_solution_1 = checksolution();
  if (num_solution_1 >= 0)
    {
      BigInt k;
      printpopulation();
      k = tobeten(population[num_solution_1].gene);
      cout << "result: gen" << generation << " cancel!!!!!!" << k << endl;
      cout << "result: N = " << N << "; 1 = " << k << "; N%k = " << N%k << "; 2 = " << N/k << endl; 
      break;
    }
    cout << generation << "  SIZE  " << population.size() << endl;
  
  if (generation == MAXGENS){
      cout << "result: gen FAILED" << endl;
      cout << "result: gen FAILED" << endl; }
  }
  // cout << "rand" << endl;
  // srand(time(NULL));
  // cout<< "0" << rand()%100/100.0 << endl;
  // cout << "1" << rand()%100/100.0 << endl;
  return 0;
}