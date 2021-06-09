# include <stdio.h>
# include <iostream>
# include <gmp.h>
# include <vector>
#include <thread>
# include <cmath>
# include <algorithm>
# include <vector>
# include <gmp.h>
# include <math.h>
# include <cassert>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <random>
# include <ctime>
# include <cstring>
# include <map>

#include <memory>

using namespace std;

# define POPSIZE 1000
# define MAXGENS 3000
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

  bool is_prime() const {
    return mpz_probab_prime_p(v, 40) > 0;
  }

  void set_next_prime(){
    mpz_nextprime(v,v);

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

struct genotype {
  vector<BigInt> gene;
  BigInt fitness;
  BigInt rfitness;
  BigInt cfitness;

  genotype() {
    fitness = 0;
    rfitness = 0.0;
    cfitness = 0.0;
  }
};

bool sort_by_fitness( const genotype & lhs, const genotype & rhs ) {
   return lhs.fitness < rhs.fitness;
}
bool compare_by_fitness( const genotype & lhs, const genotype & rhs ) {
   return lhs.fitness == rhs.fitness;
}


class common{
  public:


vector<struct genotype> population;
vector<struct genotype> newpopulation;

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
  return vec;
}

void initialize (BigInt N, int number_of_flow){
  int i;
  int j, length = 0;
   map <int, int> probability;
  vector<BigInt> Nbit = tobit(N);
  int size_chrom, sum;
  BigInt zero, first_chrom_in_ten, N_root;
   unsigned int two = 2;
  vector<BigInt> num;
  N_root = N.root(two);
  cout << N_root << endl;
  while (N_root > zero) {
    N_root /=10;
    length++;
  }
  cout << " length " << length << endl;
  
  for (int k = 2; k <= (length -2) ; k++){
    probability[k] = POPSIZE * 40 / (length - 2) / 100;
  }
  probability[length-1] = probability[length] = POPSIZE*30/100;
  cout << "probability"  << endl;
  for (const auto& [digit_size, count] : probability) {
     cout << "count " << count << "digit " << digit_size <<  endl;
     sum += count;
  }
  int need_to_add =  - sum;
  if (need_to_add != 0){
    probability[length + 1] = need_to_add;
  }
  for (const auto& [digit_size, count] : probability) {
      if (digit_size != length + 1){
          int low_bonder = pow(10, digit_size-1);
          int high_bonder = pow(10, digit_size);
        
          for (int i = 0; i < count ; i++){
          BigInt num_2;
              num_2 = rand() % (high_bonder-low_bonder) + low_bonder;
              num_2.set_next_prime();
              num.push_back(num_2);}
        }
    else{
        for (int i = 0; i < count ; i++){
            BigInt num_1;
            num_1 = rand()%((int)pow(10, (length - 2)) - 200) + 200;
            num_1.random_num();
            // cout << "num_1 : " << num_1 << endl;
            num_1.set_next_prime();
            num.push_back(num_1);}}
  }
  for ( j = 0; j < num.size(); j++ )
  {
    genotype gt;

//       // cerr << population.size() << endl;
      gt.fitness = 0;
      gt.rfitness = 0;
      gt.cfitness = 0;
      gt.gene = tobit(num[j]);
      population.push_back(gt);
    }
  cout << "population.size " << population.size() << endl;
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

// расчёт фитнесс-функции для каждой хромосомы в популяции
// void evaluate (BigInt N, int flow, int part)
// {
//   int gen;
//   int member, diff;
//   vector<BigInt> x;
//   if (population.size()%flow != 0 ){
//       cout << "population.size()/flow" << population.size() << " , " << population.size()/flow << endl;
//       diff = population.size() - flow*(population.size()/flow);
//   }
//   cout << "diff: " << diff << endl;
//   int low_b = (part - 1)*population.size()/flow;
//   cout << "low_b: " << low_b << endl;
//   int high_b = part * population.size()/flow;
//   cout << "high_b: " << high_b << endl;
//   for (member = low_b; member < high_b ; member++ )
//   {
//     BigInt myNumInInteger;
//     BigInt zero;
//     x = population[member].gene;
//     myNumInInteger = tobeten(x);
//     // cout << "myNumInInteger" << myNumInInteger << endl;
//     if (!(myNumInInteger == zero) ){
//     population[member].fitness = N% myNumInInteger;
//     }
//     else { population.erase(population.begin() + member); }
//     // cout << "population[member].fitness " << population[member].fitness << endl;
//   }
//   cout << "member" << member << endl;
//   return;
// }

void printpopulation() 
{
    cout << population.size() << endl;
    for (int i = 0; i< population.size(); i++)
    {

      cout << i << " genes value " << endl;
      cout << "in ten " << tobeten(population[i].gene) << " in bit " << endl;  
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
vector<BigInt> made_crom_with_prop(vector<BigInt> num){
  // BigInt new_num;
  // cout << tobeten(num) << endl;
  auto new_num = tobeten(num);
  new_num.set_next_prime();
  // cout << "NEW CHROM CHECK" << endl;
  // cout << tobeten(num) << endl;
  // cout << new_num << endl;
  return tobit(new_num);
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


void selector ()
{
  vector<genotype> newpop;
  // srand(time(NULL));
  int i;
  int j;
  int mem;
  double p;
  BigInt sum;

  for ( i = 0; i < population.size(); i++ )
  { 
    p = rand()%100/100.0; // [0,1]
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
        cout << " one " << one << "mem " << mem << endl; 
        // cout << "Xover начался" << endl;
        Xover ( one, mem, N);
        cout << "Xover закончился" << endl;
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
  // cout << "OLD" << endl;
  // cout << "one" << tobeten(population[one].gene) << endl;
  // cout << "two" << tobeten(population[two].gene) << endl;
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
  popmutateone.gene = made_crom_with_prop(num1);
  vector<BigInt> num2 = popmutatetwo.gene;
  auto g =  made_crom_with_prop(num2);
  popmutatetwo.gene = g;
  // cout << "обмен генами закончен" << endl;
  // if (checkpropety(tobeten(popmutateone.gene))){
  //   population.push_back(popmutatetwo);
  // }
  // else if (checkpropety(tobeten(popmutatetwo.gene))){
  //    population.push_back(popmutateone);
  // }
  // evaluate_for_one(popmutatetwo);
  population.push_back(popmutatetwo);
  // cout << "popmutateone : " << tobeten(popmutateone.gene) << endl;
  // cout << "popmutatetwo : " << tobeten(popmutatetwo.gene) << endl; 
  // << "popmutatetwo : " <<  popmutatetwo.gene.size() << endl;
  // cout << "FIRST" << tobeten(popmutateone.gene) << endl;
  population.push_back(popmutateone);
  // cout << "SECOND" << tobeten(popmutatetwo.gene) << endl;
  return;
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
  // printpopulation();
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
      // cout << i << endl;
      // cout << population.size();
      // cout << population[i].gene.size() << endl;
      vector<BigInt> copy  = population[i].gene;
      // mpz_xor(res, copy[random_number].getValue(), one.getValue());
      // cout << random_number << " check "<< copy.size() << endl;
      copy[random_number] = copy[random_number].xor_num(one);
      population[i].gene = made_crom_with_prop(copy);
      // cout << "NEW AFTER MUTATE   " << tobeten(population[i].gene) << endl;
  }
  }
  // cout << "количество мутирующих " << k << endl;
  return;
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
      maxim = N.root(two);
      while (t == zero | t == one){
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

void Run(BigInt N, int number_of_flow) {

  // TODO N
  int generation;
  initialize(N, number_of_flow);
//   printpopulation();
  evaluate(N);
  printpopulation();
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
  // return ;
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
  // может быть стоит оставить одинаковые
  population.erase(std::unique(population.begin(), population.end(), compare_by_fitness), population.end());
  // printpopulation();
  add_new_chrome(N);
  // printpopulation();
  check();
  // cout << "NEWPOPILATION SIZE: " << newpopulation.size() << endl;
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
}
};

void Run(string number, int number_of_flow) {
    BigInt N(number);
    common c;
    c.Run(N, number_of_flow);
}

int main() {
  srand(time(NULL));
  string N = "6984871";
//   //file_name - для проверок
// //   string filename = "initial_number.txt";
//   common c,d;
//   c.Run(N, filename); 
  int number_of_flow = 2;
  thread th1(Run, N, number_of_flow);
  thread th2(Run, N, number_of_flow);

  th1.join();
  th2.join();
  return 0;
}