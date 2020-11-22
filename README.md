# mathfunction
C++で&lt;cmath>ヘッダファイルで定義されていないような数学関数を入れていきたい。速度は関数によっては他の外部ライブラリで提供されているものに比べて劣るとは思うけど少しテストしたい時とかに使えるかも。一番はプロが作ったライブラリを使うこと。

・関数リスト  

```c++
double power(double x, int N);
```
標準で、std::pow()があるが指数が比較的小さな整数の時は直接N回掛けた方が速いのかなと思い導入。
```c++
double lower_incomp_gamma(double a, double x);
```
第1種不完全ガンマ関数γ(a, x)を算出する。実装は https://en.wikipedia.org/wiki/Incomplete_gamma_function にある冪級数展開の式を使用した。(故に速度が遅い)
```c++
double normalized_lower_incomp_gamma(double a, double x);
```
第1種不完全ガンマ関数をΓ(a)で正規化した値を返す。
```c++
double incomp_beta(double x, double a, double b);
```
正則化不完全ベータ関数I_x(a, b)の値を返す。実装は BOOSTの不完全ガンマ関数のページ   https://www.boost.org/doc/libs/1_50_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_beta/ibeta_function.html  
およびページ先の参考論文および下記URL  
https://www.seijo.ac.jp/pdf/faeco/kenkyu/118/118-sekimoto.pdf   
を参考にした。
```c++
double beta(double a, double b)
```
ベータ関数B(a, b)の値を返す。C++17から実装されているが、コンパイラが対応していなかったので泣く泣く実装。

・以下は統計学でよく用いられる確率関数などを実装した。  
```c++
double chisquared_pdf(double x, double deg);
```
自由度degのカイ二乗分布の確率密度関数のxにおける値を返す。

```c++
double chisquared_cdf(double x, double deg);
```
自由度degのカイ二乗分布の累積分布関数 F(x;deg)を算出する。  
```c++
double chisquared_lower_limit(double alpha, double deg, double init = 0);
```
カイ二乗分布の下側確率(下側α点)を計算する。極端な値だと収束しない可能性があるけど大まかには収束する。
initは探索の初期値でこれを変えると収束しやすくなるかもしれないが特に変える必要はない。  
```c++
double poisson_pdf(unsigned int k, double lambda);
```
平均 lambda のポアソン分布のX=kにおける確率を返す。ポアソン分布は、ある一定時間内に平均してlambda回起こる事象Aがあったとき、Aが実際に起きる回数が従う分布である。 
これは起こる確率が低い事象に対してよく当てはまり、例えばある交差点での1日の事故件数などがそれにあたる。  
```c++
double binomial_pdf(unsigned int n, unsigned int k, double p);
```
ある確率pの事象をn回試行したとき、k回起こる確率(いわゆる二項分布)を返す。表裏が確率0.5で出るコインを10回投げた時、表が5回出る確率、など。

```c++
double binomial_cdf(unsigned int n, unsigned int k, double p);
```
二項分布の累積和を返す。ある確率pの事象Aをn回試した時、k回以下Aが起こる確率を表す。例えば同様に確からしいコインを100回投げて50回以下表が出る確率、など。

```c++
double normal_pdf(double x, double mu = 0, double sigma = 1);
```
平均mu, 標準偏差sigmaにおける正規分布の確率密度関数 f(x) を返す。
```c++
double normal_cdf(double x, double mu = 0, double sigma = 1);
```
平均mu, 標準偏差sigmaにおける正規分布の累積分布関数 F(x) を返す。
