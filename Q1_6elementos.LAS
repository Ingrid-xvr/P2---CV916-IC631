//Questao 1 - P2 CV916 

// Coordenadas dos nos
Int XY = [	0	0
		1/6	0
		2/6     0
		3/6     0
		4/6     0
	 	5/6     0
    		1   0]

// Criacao dos elementos
Int Beam_EL = [	1	2	
			2     3
			3	4
			4	5
			5	6
			6	7]	
			  
// Definicao dos graus de liberdade para cada no
Int EQ = [	0 0 1
		0 2 3
		0 4 5
		0 6 7
		0 8 9
		0 10 11
  		0 0  12]

// Locacao dos elementos na matriz
Int Beam_LM = Beam_Make_LM(Beam_EL, EQ) 

// Propriedades das matrizes [A As I E Nu Dn]
Dbl P1 = [1.0 0.0 1.0 1.0 0.0 0.0]
Dbl Prop = [1.0 0.0 1.0 1.0 0.0 0.0] 

// Definicao da matriz de rigidez local
Dbl KE1 = Beam_Make_K(1,XY,Beam_EL,P1,"Global")
Dbl KE2 = Beam_Make_K(2,XY,Beam_EL,P1,"Global") 
Dbl KE3 = Beam_Make_K(3,XY,Beam_EL,P1,"Global")  
Dbl KE4 = Beam_Make_K(4,XY,Beam_EL,P1,"Global")
Dbl KE5 = Beam_Make_K(5,XY,Beam_EL,P1,"Global") 
Dbl KE6 = Beam_Make_K(6,XY,Beam_EL,P1,"Global") 

// Montagem da matriz de rigidez global
Dbl KE = Zero (12) 
Assemble (KE, KE1, Beam_LM, 1)
Assemble (KE, KE2, Beam_LM, 2)
Assemble (KE, KE3, Beam_LM, 3)
Assemble (KE, KE4, Beam_LM, 4)
Assemble (KE, KE5, Beam_LM, 5)
Assemble (KE, KE6, Beam_LM, 6)

// definicao da matriz de rigidez geométrica local
//KG = Beam_Make_KG(N ,e , XY, Beam_EL, Prop, Ref, Type)
Dbl KC11 = Beam_Make_KG(1,1,XY,Beam_EL,Prop,"Global","C")
Dbl KC22 = Beam_Make_KG(1,2,XY,Beam_EL,Prop,"Global","C")
Dbl KC33 = Beam_Make_KG(1,3,XY,Beam_EL,Prop,"Global","C")
Dbl KC44 = Beam_Make_KG(1,4,XY,Beam_EL,Prop,"Global","C")
Dbl KC55 = Beam_Make_KG(1,5,XY,Beam_EL,Prop,"Global","C")
Dbl KC66 = Beam_Make_KG(1,6,XY,Beam_EL,Prop,"Global","C")

// montagem da matriz de rigidez geométrica global
Dbl KG = Zero (12)
Assemble (KG, KC11, Beam_LM, 1)
Assemble (KG, KC22, Beam_LM, 2)
Assemble (KG, KC33, Beam_LM, 3)
Assemble (KG, KC44, Beam_LM, 4)
Assemble (KG, KC55, Beam_LM, 5)
Assemble (KG, KC66, Beam_LM, 6)

Dbl {Modos,lambda} = Jacobi(KG,KE)

End
