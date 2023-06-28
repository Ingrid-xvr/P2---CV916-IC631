//Questao 1 - P2 CV916

// Coordenadas dos nos
Int XY = [	0	0
		1/3	0
		2/3     0
    		1   0]

// Criacao dos elementos
Int Beam_EL = [	1	2	
			2     3
			3	4]	
			  
// Definicao dos graus de liberdade para cada no
Int EQ = [	0 0 1
		0 2 3
		0 4 5
		0 0 6]

// Locacao dos elementos na matriz
Int Beam_LM = Beam_Make_LM(Beam_EL, EQ) 

// Propriedades das matrizes [A As I E Nu Dn]
Dbl P1 = [1.0 0.0 1.0 1.0 0.0 0.0]
Dbl P2 = [1.0 0.0 1.0 1.0 0.0 0.0]
Dbl P3 = [1.0 0.0 1.0 1.0 0.0 0.0]
Dbl Prop = [1.0 0.0 1.0 1.0 0.0 0.0] 

// Definicao da matriz de rigidez local
Dbl KE1 = Beam_Make_K(1,XY,Beam_EL,P1,"Global")
Dbl KE2 = Beam_Make_K(2,XY,Beam_EL,P2,"Global") 
Dbl KE3 = Beam_Make_K(3,XY,Beam_EL,P3,"Global")  

// Montagem da matriz de rigidez global
Dbl KE = Zero (6) 
Assemble (KE, KE1, Beam_LM, 1)
Assemble (KE, KE2, Beam_LM, 2)
Assemble (KE, KE3, Beam_LM, 3)


// definicao da matriz de rigidez geométrica local
//KG = Beam_Make_KG(N ,e , XY, Beam_EL, Prop, Ref, Type)
Dbl KC11 = Beam_Make_KG(1,1,XY,Beam_EL,Prop,"Global","C")
Dbl KC22 = Beam_Make_KG(1,2,XY,Beam_EL,Prop,"Global","C")
Dbl KC33 = Beam_Make_KG(1,3,XY,Beam_EL,Prop,"Global","C")

// montagem da matriz de rigidez geométrica global
Dbl KG = Zero (6)
Assemble (KG, KC11, Beam_LM, 1)
Assemble (KG, KC22, Beam_LM, 2)
Assemble (KG, KC33, Beam_LM, 3)


End
