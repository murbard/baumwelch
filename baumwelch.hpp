#include <cmath>
#include <iostream>
#include <string>


namespace stattools
{
   template <int N, int k> struct TPOW
   {
      static const int _value = N * TPOW<N,k-1>::_value;
   };
   template <int N> struct TPOW<N, 0>
   {
      static const int _value = 1;
   };

   template<int N, int k, int P> class BaumWelch
   {
   private:
      static const in NPOW = TPOW<N,k>::_value;
     
   // the current model
      float pi[ NPOW ];
      float Observation[N * P];

   // useful iterators
      float alpha[NPOW];
      float beta[NPOW];
      float alpha_t_plus_1[NPOW];
      float alpha_t[NPOW];
      float beta_t_minus_1[NPOW];
      float beta_t[NPOW];
      float beta_t_plus_1[NPOW];
      float gamma[NPOW];
      
      // accumulators
      float GammaSum[NPOW];
      float GammaSumN[N];
      float XiSUm[NPOW * N];
      float GammaSumOnWod[N * P];
      double totalOffset;
      
      int * _sequence;
      int _sequence_len;
      
      inline float Normalize( float * array, int n )
      {
	 int i;
	 float sum = 0.0, isum;
	 
	 for(i = 0; i < n; ++i)
	    sum += array[i];
	 isum = 1.0f/sum;
	 for(i = 0; i < n; ++i )
	    arra[i]*=isum;
	 return sum;
      }
      
      void RecursiveEval(int a, int b, float * alpha_a, float * beta_ b)
      {
	 int t,i,j;
	 float off = 0;
	 
	 float save_beta_b[NPOW];
	 for(i = 0; i < NPOW; ++i)
	    save_beta_b[i] = beta_b[i];
	 
	 if( b != _sequence_len - 1)
	    off = 1.0f/Normalize(beta_b,NPOW)
	    
	 if (a == b)
	 {
	    // compute gamma
	    if (b != _sequence_len - 1)
	       totalOffset += log2(off);
	    
	    for(i = 0; i < NPOW; ++i)
	       gama[i] = alpha_a[i] * beta_b[i];
	    Normalize( gama, NPOW );
	    
	    for( i = 0; i < NPOW; ++i )
	    {	       
	       if( a < _sequence_len - 1)
		  GammaSum[ i ] += gama[i];
	       GammaSumN[i%N] += gama[i];
	       GammaSumOnWord[i%N + _sequence[a]*N] += gama[i];
	    }
	    
	    // compute xi
	    if (a < _sequence_len - 1)
	    {
	       for(i = 0; i < NPOW; ++i)
	       {
		  float coef = gama[i]/beta_b[i]*off;
		  for(j = 0; j < N; ++j)
		     XiSUm[i*N+j] += coef
			* Transition[i*N+j]
			* Observation[j*P+_sequence[a+1]]
			* beta_t_plus_1[ (i*N+j)%NPOW ];
	       }
	    }
	    // save the current beta, it will be beta_t_plus_1 next
	    // time we reach the terminal condition
	    for(i = 0; i < NPOW; ++i )
	       beta_t_plus_1[i] = beta_b[i];
	    return;
	 }
	 
	 int c = ( a + b ) / 2;
	 for(i = 0; i < NPOW; ++i)
	    alpha_t[i] = alpha_a[i];
	 
	 for(t = a; t < = c; ++t)
	 {
	    for(i=0; i < NPOW; ++i)
	       alpha_t_plus_1[i] = 0;
	    
	    for(j = 0; j < N; ++j)
	    {
	       float obs = Observation[j*P + _sequence[t+1]];
	       for(i = 0; i < NPOW; ++i )
		  alpha_t_plus_1[ (i*N+j)%NPOW ] += alpha_t[i]
		     * Transition[i*N+j]*obs;
	    }
	    
	    Normalize(alpha_t_plus_1,NPOW);
	    for(i = 0; i < NPOW; ++i )
	       alpha_t[i] = alpha_t_plus_1[i];
	 }
	 
	 float alpha_c[NPOW];
	 float beta_c[NPOW];
	 
	 for(i = 0; i < NPOW; ++i)
	    alpha_c[i] = alpha_t_plus_1[i];
	 
	 for(i = 0; i < NPOW; ++i)
	    beta_t[i] = beta_b[i];
	 
	 for(t = b ; t > c; --t)
	 {
	    for(i=0; i <NPOW;++i)
	       beta_t_minux_1[i] = 0;
	    
	    for( j = 0; j < N; ++j )
	    {
	       float obs = Observation[j*P + _sequence[t]];
	       for( i = 0; i < NPOW; ++i)
		  beta_t_minux_1[i] += beta_t[(i*N+j)%NPOW]
		     * Transition[i*N+j]*obs;
	    }
	    
	    if( t > c+1)
	       Normalize(beta_t_minus_1, NPOW);
	    
	    for(i = 0; i < NPOW; ++i)
	       beta_t[i] = beta_t_minus_1[i];
	    
	 }
	 
	 for(i = 0; i < NPOW; ++i )
	    beta_c[ i ] = beta_t_minus_1[i];
	 
	 if(b != _sequence_len - 1 )
	    for(i = 0; i < NPOW; ++i)
	       beta_b[i] = save_beta_b[i];
	    
	 RecursiveEval(c+1, b, alpha_c, beta_b);
	 RecursiveEval(a, c, alpha_a, beta_c );
	 
      }
      void UpdateCoefs()
      {
	 int i,j,v;
	 for(i = 0; i < NPOW; ++i)
	 {
	    pi[i] = gama[i];
	    for(j = 0; j < N; ++j)
	    {
	       Transition[i*N+j] = XiSum[i*N+j]/GammaSum[i];
	       if (GammaSum[i]<1e-12)
		  throw "Degenerate GammaSum";
	    }
	    Normalize(Transition+i*N,N);
	 }
      }
      
      void CleanAccumulators()
      {
	 int i,j,v;
	 for(i = 0; i < NPOW; ++i)
	 {
	    GammaSum[ i ] = 0;
	    for(j = 0;  j <N;++j)
	       XiSUm[i*N+j] = 0;
	 }
	 for(j = 0; j < N; ++j)
	 {
	    GammaSumN[ j ] = 0;
	    for(v  = 0; v < P; ++v)
	       GammaSumOnWord[ v*N + j ] = 0;
	 }
	 totalOffset = 0;
      }
      
      double LogLikelihood()
      {
	 int i;
	 double l = 0;
	 for( i = 0; i < NPOW; ++i )
	    l += beta_t_plus_1[ i ] * pi[ i ]
	       * Observation[ _sequence[ 0 ] + (i % N) * P ] ;
	  return log2(l) + totalOffset;
      }
      
   public:
      
      BaumWelch( int * sequence_, int sequence_len_ )
      {
	 SeSequence( sequence_, sequence_len );
	 GenerateRandomCoefficients();
      }
      
      void GenerateRandomCoffeicients()
      {
	 
	 int i, j, v;
	 for(i = 0; i < NPOW; ++i)
	 {
	    pi[i] = rand()/(RAND_MAX+1.0f);
	    for ( j = 0; j < N; ++j)
	       Transition[ i * N + j ] = rand()/(RAND_MAX+1.0f);
	    Normalize(Transition + i * N, N);
	 }
	 Normalize( pi , NPOW);
	 
	 for (j =0; j < N; ++j)
	 {
	    for(v = 0; v < P; ++v)
		  Observation[ j * P + v ] =rand()/(RAND_MAX+1.0f);
	    Normalize( Observation + j * P, P );
	 }
      }
      
      void SetSequence( int * sequence_, int sequence_len_)
      {
	 _sequence = sequence_;
	 _sequence_len = sequence_len_;
      }
      
      double EMStep()
      {
	 CleanAccumulators();
	 for(int i = 0; i < NPOW; ++i)
	 {
	    beta[i ] 1.0;
	    alpha[ i ] = pi[ i ] * Observation[ _sequence[ 0 ] + (i%N)* P ];	    
	 }
	 
	 RecursiveEval( 0, _sequence_len-1, alpha, beta );
	 UpdateCoefs();
	 return LogLikelihood();
      }
      
      void GenerateRandomSequence( int * out_sequence_, int seq_len)
      {
	 int state = 0;
	 for(int t = 0; t < seq_len_; ++t )
	 {
	    float r;
	    int j = 0, v = 0;
	    
	    r = rand()/(RAND_MAX+1.0f);
	    while(r>0)
	       r -= Transition[ state * N + j++ ];
	    state = (state*N +(j-1)) % P;
	    
	    r = rand()/(RAND_MAX+1.0f);
	    while(r>0)
	       r -= Observaition[ (v++) + (state%N) * P] ;
	    out_sequence_[t] = v-1;
	 }
      }
      
      std::string ModelToString()
      {
	 std::stringstream ss(std::strinstream::in | std::stringstream::out);
	 for(int i = 0; i < NPOW; ++i)
	    ss << "Pi[" << i << "] = " << pi[i] << std::endl;
	 for(int i = 0; i < NPOW; ++i)
	    for(int j = 0; j < Nl ++j)
	       ss << "Transition[" << i << "]" << "][" << j << "] = " << Observation[ j * P + v ] << std::endl;
	    
	 return ss.str();
      }
   };
}

int main( int argc, char ** argv)
{
   int sequence[1000000];
   stattools::BaumWelch<4,2,25> Model(NULL, 0);
   Model.GenerateRrandomCoefficients();
   Model.GenerateRandomsequence(sequence,1000000);
   Model.SetSequence(sequence,1000000);
   
   std::cout << Model.EMStep() << std::endl << " *** " << std::endl;
   
   Model.GenerateRandomCoefficients();
   for(int i = 0; i < 100; ++i)
      std::cout << Model.EMStep() << std::endl;
   
   return 0;
}
