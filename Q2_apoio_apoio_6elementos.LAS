//Questao 2

// Coordenadas dos nos
Int XY = [	0	0
            0     1/6
            0     2/6
            0     3/6
            0     4/6
		0	5/6
	      0	1]

// Criacao dos elementos
Int Beam_EL = [	1	2 
                  2     3
                  3     4
                  4     5
			5	6
		      6	7]	
			  
// Definicao dos graus de liberdade para cada no
Int EQ = [	0 0 1
		2 0 3
            4 0 5
            6 0 7
            8 0 9
		10 0 11
		0 0 12]

// Locacao dos elementos na matriz
Int Beam_LM = Beam_Make_LM(Beam_EL, EQ) 

// Propriedades das matrizes [A As I E Nu Dn]
Dbl P1 = [0.053 0.0 5.359*10^(-4) 2.3*10^(10) 0.0 0.0]

// Definicao da matriz de rigidez local
Dbl KE1 = Beam_Make_K(1,XY,Beam_EL,P1,"Global")
Dbl KE2 = Beam_Make_K(2,XY,Beam_EL,P1,"Global")
Dbl KE3 = Beam_Make_K(3,XY,Beam_EL,P1,"Global")
Dbl KE4 = Beam_Make_K(4,XY,Beam_EL,P1,"Global")
Dbl KE5 = Beam_Make_K(5,XY,Beam_EL,P1,"Global")

// Montagem da matriz de rigidez global
Dbl K = Zero (12) 
Assemble (K, KE1, Beam_LM, 1)
Assemble (K, KE2, Beam_LM, 2)
Assemble (K, KE3, Beam_LM, 3)
Assemble (K, KE4, Beam_LM, 4)
Assemble (K, KE5, Beam_LM, 5)
Assemble (K, KE5, Beam_LM, 6)

//Matrizes de incidencia
Dbl BB1 = [0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0]

Dbl BB2 = [0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0]

Dbl BB3 = [0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0]

Dbl BB4 = [0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0]

Dbl BB5 = [0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0]

Dbl BB6 = [0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0]


Dbl E = 2.3*10^(10)
Dbl I = 5.359*10^(-4)
Dbl L = 1/6
Dbl aux = E*I/L^3

//Matriz de rigidez elementar
Dbl K_local = aux*[12  6*L  -12  6*L                
			     6*L 4*L^2 -6*L 2*L^2               
			     -12 -6*L 12   -6*L                
			     6*L 2*L^2 -6*L 4*L^2]

//Matriz de rigidez elastica global
Dbl K_E = t(BB1)*K_local*BB1+t(BB2)*K_local*BB2+t(BB3)*K_local*BB3+t(BB4)*K_local*BB4+t(BB5)*K_local*BB5+t(BB6)*K_local*BB6

Dbl rho = 2500
Dbl A = 0.053

//Matriz de massa elementar
Dbl M_local = rho*A*L/420*[156	22*L  54  -13*L                
			     	 22*L 4*L^2 13*L -3*L^2               
			       54	13*L  156   -22*L                
			      -13*L -3*L^2 -22*L 4*L^2]

//Matriz de rigidez elastica global
Dbl M = t(BB1)*M_local*BB1+t(BB2)*M_local*BB2+t(BB3)*M_local*BB3+t(BB4)*M_local*BB4+t(BB5)*M_local*BB5+t(BB6)*M_local*BB6
		
// Analise modal
Dbl {Phi,W2} = Jacobi(M,K_E)
Dbl W = Sqrt(W2)  //rad/s
Dbl Periods= W@^(-1) @* 2*3.141592654 //s
Dbl Freq = W@/(2*3.14159265) //Hz

End
