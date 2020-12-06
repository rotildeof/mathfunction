# mathfunction
C++で&lt;cmath>ヘッダファイルで定義されていないような数学関数を入れていきたい。速度は関数によっては他の外部ライブラリで提供されているものに比べて劣るとは思うけど少しテストしたい時とかに使えるかも。一番はプロが作ったライブラリを使うこと。
名前空間はmathfunc。

数値計算用関数
- 

---
```c++
double power(double x, int N);
```
標準で、std::pow()があるが指数が比較的小さな整数の時は直接N回掛けた方が速いのかなと思い導入。

---
```c++
double differential(std::function<double(double)> func, double x, double h = 0.001);
```

```c++
double differential2(std::function<double(double)> func, double x, double h = 0.001);
```
(differencial) : 連続関数f(x)のxにおける微分係数f'(x)を求める。hは刻み幅でデフォルトで0.001。この幅によっては精度が変わる。
(differencial2) : 連続関数f(x)のxにおける二階微分f''(x)を求める。

```c++
// 例 (#includeなど適宜補完してください)
auto f = [](double x){return std::sin(x);}; // f(x) = sin(x)
double dif1 = mathfunc::differential(f, 0) // f'(0)
double dif2 = mathfunc::differential2(f, 0) // f''(0)
std::cout << dif1 << std::endl;
std::cout << dif2 << std::endl;
  
// 出力
// 1.0000000
// 0.0000000
```

---

```c++
double error_propagation(std::function<double(double*)> func, double* x, double* x_e, const int num_arg, double h = 0.001);
```

誤差伝搬用の関数。例えばある測定値Xと測定値Yの(統計)誤差がそれぞれa, bのとき、Z=f(X, Y)としたときのZの(統計)誤差を求める。

```c++
// 例 z = (x + y) / 2 の場合
auto f = [](double *x){return (x[0] + x[1]) / 2;} // z = (x + y) / 2;
double x[2] = {2, 3}; // 測定値がx = 2, y = 3 だったとする
double x_e[2] = {0.1, 0.1}; // x, yの測定誤差がそれぞれ 0.1, 0.1 だったとする。
double result = mathfunc::error_propagation(f, x, x_e, 2); // 第3引数はZ = f(X, Y,...) の fの引数の数。今の場合は x と y の2つ。
std::cout << result << std::endl; // z の(統計)誤差の結果を算出する。

// 出力
// 0.070710678
```

数値計算アルゴリズム
-

---
```c++
double newton_method(std::function<double(double)> func, double init, double epsilon = 1e-12);
```

方程式 f(x)=0 の近似解をニュートン法を用いて求める。initで与える値によっては収束せず、解の近くの値を与える必要がある。

```c++
// 例
auto f = [](double x){return x * x - 2;}; // f(x) = x^2 - 2 
double result = mathfunc::newton_method(f, 1) // 初期値 1 で f(x) = 0 を解く。
std::cout << result << std::endl;

// 出力 : ルート2に近い値が得られる。
// 1.4142136
```

---
```c++
double find_extremum_x(std::function<double(double)> func, double init, double epsilon = 1e-12);
```
関数 y=f(x)が極値をとるxの近似値を求める。initの近くの極小値または極大値を求める。initによっては収束しない。

```c++
// 例
auto f = [](double x){return std::sin(x);}; // f(x) = sin(x)
double result = mathfunc::newton_method(f, 1) // 初期値 1 付近で sin(x) が極値をとるxを探索する。
std::cout << result << std::endl;

// 出力 : pi/2 に近い値が得られる。
// 1.5707963

```

---
```c++
int64_t gcd(int64_t a, int64_t b);
int64_t lcm(int64_t a, int64_t b);
```

(gcd) : 整数 a, b の最大公約数を求める。
(lcm) : 整数 a, b の最大公倍数を求める。

特殊関数
-

---
```c++
double lower_incomp_gamma(double a, double x);
```

第1種不完全ガンマ関数γ(a, x)を算出する。実装は https://en.wikipedia.org/wiki/Incomplete_gamma_function にある冪級数展開の式を使用した。(故に速度が遅い)

---
```c++
double normalized_lower_incomp_gamma(double a, double x);
```
第1種不完全ガンマ関数をΓ(a)で正規化した値を返す。

---
```c++
double incomp_beta(double x, double a, double b);
```
正則化不完全ベータ関数I_x(a, b)の値を返す。実装は BOOSTの不完全ガンマ関数のページ   https://www.boost.org/doc/libs/1_50_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_beta/ibeta_function.html  
およびページ先の参考論文および下記URL  
https://www.seijo.ac.jp/pdf/faeco/kenkyu/118/118-sekimoto.pdf   
を参考にした。(この連分数の導出法知ってる人いたらどなたか教えて欲しい・・・。)

---
```c++
double beta(double a, double b)
```
ベータ関数B(a, b)の値を返す。C++17から実装されているが、コンパイラが対応していなかったので泣く泣く実装。

統計学で使う関数
-

---
```c++
double chisquared_pdf(double x, double deg);
```
自由度degのカイ二乗分布の確率密度関数のxにおける値を返す。

---
```c++
double chisquared_cdf(double x, double deg);
```
自由度degのカイ二乗分布の累積分布関数 F(x;deg)を算出する。  

---
```c++
double chisquared_lower_limit(double alpha, double deg, double init = 0);
```
カイ二乗分布の下側確率(下側α点)を計算する。極端な値だと収束しない可能性があるけど大まかには収束する。
initは探索の初期値でこれを変えると収束しやすくなるかもしれないが特に変える必要はない。  

---
```c++
double geometric_pdf(int k, double p);
```
確率pのベルヌーイ試行を繰り返した時、X=kで初めて成功する確率を求める(Xが従う分布を幾何分布という)。例えば、サイコロを振ってk=5回目に初めて6の目が出る確率、など。

---
```c++
double geometric_cdf(int k, double p);
```
幾何分布の累積分布関数を求める。

---
```c++
double poisson_pdf(unsigned int k, double lambda);
```
平均 lambda のポアソン分布のX=kにおける確率を返す。ポアソン分布は、ある一定時間内に平均してlambda回起こる事象Aがあったとき、Aが実際に起きる回数が従う分布である。 
これは起こる確率が低い事象に対してよく当てはまり、例えばある交差点での1日の事故件数などがそれにあたる。  

---
```c++
double binomial_pdf(unsigned int n, unsigned int k, double p);
```
ある確率pの事象をn回試行したとき、k回起こる確率(いわゆる二項分布)を返す。表裏が確率0.5で出るコインを10回投げた時、表が5回出る確率、など。

---
```c++
double binomial_cdf(unsigned int n, unsigned int k, double p);
```
二項分布の累積和を返す。ある確率pの事象Aをn回試した時、k回以下Aが起こる確率を表す。例えば同様に確からしいコインを100回投げて50回以下表が出る確率、など。
(実装は実直に足すのではなく、不完全ベータ関数で算出している。)

---
```c++
double normal_pdf(double x, double mu = 0, double sigma = 1);
```
平均mu, 標準偏差sigmaにおける正規分布の確率密度関数 f(x) を返す。

---
```c++
double normal_cdf(double x, double mu = 0, double sigma = 1);
```
平均mu, 標準偏差sigmaにおける正規分布の累積分布関数 F(x) を返す。

---
