function [ind v] = simplex_res(A,b,c,m,n,x,indB,Binv)

	fprintf('\nSimplex: Fase 2 -  Método Revisado\n')

	# Condição para entrar no while

   contador = 1;
   
	# Loop principal do programa
	while contador == 1  || min(c_barra) < 0
		

    
		# Primeira iteração
		fprintf('\n')	
		fprintf('***************** ITERANDO %d *****************\n', contador-1)
		fprintf('\n')	
		
		# Printa índices e valores das variáveis básicas
		fprintf('Índices e valores das variáveis básicas:\n')
		for i = 1: size(indB,1)
			fprintf('%d   %f\n',indB(i),x(indB)(i))
		endfor
		fprintf('\n')
    
    # Definindo os vetores associados às variáveis básicas
		# Matriz básica B 
		B = A(:,indB); 
    
    fprintf('Matriz inversa da base:\n')
    disp(Binv)
    fprintf('\n')
    
		# Variáveis básicas associadas à solução básica inicial
		x_b = x(indB);

		# Vetor de custo associado às variáveis básicas
		c_b = c(indB);

		# Cálculo do vetor p = (c_B)^T * B^-1
    
    # variáveis não básicas
    var_nb = setdiff(1:length(x), indB);
    
		fprintf('Vetor p:\n')
		p = c_b' * Binv;
    
		for i = 1: size(p,2)
			fprintf("   %f\n", p(i))  
		endfor
		fprintf('\n')

		# Inicializa o vetor c_barra com zeros
		c_barra = zeros(1, n);

		# Faz o cálculo de vetor de custos c_j_barra = c_j - (p * A_j); em que p = (c_B)^T * B^-1
		for j = var_nb
      # Cálculo do custo reduzido para as variáveis não básicas
      c_barra(j) = c(j) - p *(A(:,j)); 
		endfor


		# Printa os custos reduzidos calculados para as variáveis não básicas, uma vez que para as variáveis básicas este custo é 0 
		fprintf('Custos Reduzidos:\n')

		# Vetor com as variáveis não básicas
    
    # faz a diferença de conjuntos para determinar os índices das variáveis não básicas
		for i = var_nb;
			fprintf('  %d   %f\n', i, c_barra(i))  
		endfor
		fprintf('\n')
    
    # Condições sobre o c_barra
		# Critério de parada -> Todos os elementos não negativos: solução ótima encontrada e o algoritmo termina
		if min(c_barra) >= 0
			ind = 0;
			v = x;
			break  
		endif
		
    # Caso contrário, se tiver algum índice negativo, ele será o índice j da direção básica viável
		j = find(c_barra < 0, 1);

		# Construção de u (1 x n), o vetor direção (j-ésima direção básica)
		u = zeros(1, n)';
    
		# Calcula u = B^-1 * A_j para as variáveis básicas
		u(indB) = Binv * A(:,j);
    
    # Condições sobre o vetor u
		# Critério de parada -> Todos os valores são <= 0: Problema inviável, custo ótimo -inf e o algoritmo termina
		if max(u(indB)) <= 0	
			ind = -1;
			v = u;
			break  
		endif
		
		# Caso contrário, ou seja, u_j > 0 para algum j, atribuimos valor 1 à j-ésima posição
		u(j) = 1;
		
		# Início do cálculo do theta_estrela, inicializando com zeros		
		theta = zeros(1, size(indB)); 
    
    
    # Calcula theta, que é o maior valor que vamos percorrer na j-ésima direção básica usando o -d = u como definição
    # Faz o cálculo do theta(i) = (x_B)(i) / (u_B)(i)
    u_basico = u(indB);
    for i = 1:size(indB, 1)
			if u_basico(i) > 0	
			  theta(i) = x_b(i)/u_basico(i); 
			endif
		endfor
    
    
    # Selecionamos o theta_minimo = min(theta), com theta != 0 e o índice associado ao valor de theta_minimo 
    theta_minimo = min(theta(theta != 0));
    indice_theta_minimo = find(theta == theta_minimo);
    
    # Printa as informações de quem entra na base, direção, theta mínimo e quem sai da base
		fprintf('Entra na base: %d\n', j)
		fprintf('\n')

		fprintf('Direção:\n')
		for i = 1: size(indB,1)
			fprintf('  %d   %f\n',indB(i),u(indB)(i))  
		endfor
		fprintf('\n')

		fprintf('Theta* = %f\n', theta_minimo)
		fprintf('\n')

		fprintf('Sai da base: %d\n', indB(indice_theta_minimo)) 
		fprintf('\n')
     
     
    # Atribui theta_estrela à posição j de x
		x(j) = theta_minimo;

	  # Faz o cálculo da nova solução viável básica para as variáveis básicas
		x(indB) = x_b - theta_minimo*u(indB);

		# T é a matriz aumentada [inv(B) , vetor u]
		T = [Binv u(indB)];

    
	  # Índices básicos: j enrra na base e indice_theta_minimo sai da base (l)
		indB(indice_theta_minimo) = j;
     
    # rotina do pivotamento da matriz T
    T(indice_theta_minimo, :) = T(indice_theta_minimo, :)/...
    T(indice_theta_minimo, size(T, 2));

    for i = 1:m
      if i != indice_theta_minimo
        T(i, :) = - T(i, size(T, 2)) * T(indice_theta_minimo, :) + T(i, :);
      endif
    endfor
    
		# Atualiza B^-1
		Binv = T(:,1:m);
	
		# Aumenta o contador
		contador = contador + 1;
		
	endwhile
	
# se tiver solução (ind = 0), retorna o custo ótimo, os índices e valores das variáveis básicas	
  if ind == 0
	  fprintf('Solução ótima encontrada com custo: %d\n', dot(c,x))
	  fprintf('\n')
	  fprintf('Índices e solução ótima:\n')
	  for i = 1: size(x,1)
		  fprintf('%d   %f\n', i, x(i))  
	  endfor
	  fprintf('\n\n\n ****** Fim de programa - Retornos da função: \n\n')

  # Se não tiver solução, retorna o custo ótimo (- infinito), índices e valores do vetor direção
  elseif ind == -1
	  fprintf('Problema ilimitado. Solução ótima não encontrada. Custo ótimo: -infinito\n')
	  fprintf('\n')
	  fprintf('Índices e vetor Direção:\n')
	  for i = 1: size(u,1)
		  fprintf('%d   %f\n', i, u(i))  
	  endfor
	  fprintf('\n\n\n ****** Fim de programa - Retornos da função: \n\n')
  endif
  
endfunction









