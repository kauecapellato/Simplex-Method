function [ind v] = simplex_tab(indB,tableau)
	
	fprintf('\n*********************** Simplex: Fase 2 -  Método Tableau ***********************\n\n')

	# Condição para entrar no while
  contador = 1;

	# while loop para a realização dos cálculos
	while contador == 1  || min(c_barra) < 0
		
     
		# Primeira iteração
		fprintf('\n')	
		fprintf('********************************* ITERANDO %d *********************************\n', contador-1)
		fprintf('\n')	
		
 
    fprintf('\n')
		# A entrada do tableau deverá ser feita por linhas, ou seja, uma linha inteira de cada vez,
    # assim este programa funcionará corretamente sem problemas de formato.
    
		# A SVB está na primeira coluna, segunda linha em diante 
		x = tableau([2:end],1);
    
    fprintf('\n')
    
		# o vetor de custos reduzidos está na primeira linha, segunda coluna em diante
		c_barra = tableau(1,[2:end]);
     
    m = size(x, 1);
    n = size(c_barra, 2);
    
    
    fprintf('Tableau:\n')
		disp(tableau)
		fprintf('\n')
    
		# Verificando o vetor de cutos reduzidos, se for >= 0 o algoritmo acaba
		if min(c_barra) >= 0
			ind = 0;
       
      # monta o vetor solução x, que começou com m variáveis e vai para n variáveis 
			solucao_x = zeros(1, size(c_barra, 2));
			solucao_x(indB) = x;
      v = solucao_x';
			break  
		endif
		
		# se tiver algum índice negativo, ele será o índice j da direção básica viável
		j = find(c_barra < 0 , 1);
	  
    
		# Construção de u (1 x n), o vetor direção (j-ésima direção básica)
		u = tableau(2:end,j+1);
    
		# Condições sobre o vetor u
		# Critério de parada -> Todos os valores são <= 0: Problema inviável, custo ótimo -inf e o algoritmo termina
		if max(u) <= 0	
			ind = -1;
      
			# monta o vetor u, que começou com m variáveis e vai para n variáveis 
			solucao_u = zeros(1, size(c_barra, 2));
			solucao_u(indB) = u;
			v = solucao_u'; 
			break  
		endif
		
		# Calcula theta, que é o maior valor que vamos percorrer na j-ésima direção básica usando o -d = u como definição
    # Faz o cálculo do theta(i) = (x_B)(i) / (u_B)(i)
    
		theta = zeros(1, m);
		for i = 1:m
			if u(i) > 0	
			  theta(i) = x(i)/u(i);
			endif
		endfor
    
    
		# Selecionamos o theta_minimo = min(theta), com theta != 0 e o primeiro índice associado ao valor de theta_minimo  
    theta_minimo = min(theta(theta != 0));
    indice_theta_minimo = find(theta == theta_minimo)(1); 
    
		# Printa as informações de quem entra na base, direção, theta mínimo e quem sai da base
		fprintf('Entra na base: %d\n', j)
		fprintf('\n')

		fprintf('Theta* = %f\n', theta_minimo)
		fprintf('\n')

		fprintf('Sai da base: %d\n', indB(indice_theta_minimo)) 
		fprintf('\n')

		# Pivoteamento - funciona de maneira semelhante ao método revisado
    
    tableau(indice_theta_minimo+1, :) = tableau(indice_theta_minimo+1, :)/tableau(indice_theta_minimo+1, j+1);
   
    for i = 1:m+1
      if i != indice_theta_minimo+1
        tableau(i, :) = - tableau(i, j+1) * tableau(indice_theta_minimo+1, :) + tableau(i, :);
      endif
    endfor
    

		# adiciona o valor j nas variáveis básicas (entrou na base) no lugar do valor indice_theta_minimo
		indB(indice_theta_minimo) = j;

		# aumenta o contador
		contador = contador + 1;
		
	endwhile
	
  
  # se tiver solução (ind = 0), retorna o custo ótimo, os índices e valores das variáveis básicas	
	if ind == 0
		fprintf('Solução ótima encontrada com custo: %d\n', -tableau(1,1))
		fprintf('\n')
		fprintf('Índices e solução viável básica ótima:\n')
		for i = 1:size(solucao_x,2)
			fprintf('%d   %f\n', i, solucao_x(i))
		endfor
		fprintf('\n\n\n ****** Fim de programa - Retornos da função: \n\n')
  
  # Se não tiver solução, retorna o custo ótimo (- infinito), índices e valores do vetor direção
	elseif ind == -1
		fprintf('Problema ilimitado. Solução ótima não encontrada. Custo ótimo: -infinito\n')
		fprintf('\n')
		fprintf('Índices e vetor Direção:\n')
		for i = 1:size(solucao_u,2)
			fprintf('%d   %f\n', i, solucao_u(i))
		endfor
		fprintf('\n\n\n ****** Fim de programa - Retornos da função: \n\n')
	endif
	
endfunction










