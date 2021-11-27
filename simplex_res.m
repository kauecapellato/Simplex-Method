function [ind x v] = simplex_res(A,b,c,m,n)

	fprintf('\n*********************** Simplex: Fase 1 -  Método Revisado ************************\n\n') 

	# Condição para entrar no while

  # faz a regra para ter b >= 0
  ind_negativo_b = find(b < 0);
  
  # multiplica as linhas de A e b [A|b] por -1 nos índices ind_negativo_b para garantir que b >= 0
  
  b(ind_negativo_b) = -b(ind_negativo_b);
  A(ind_negativo_b, 1:end) = -A(ind_negativo_b, 1:end);
  
  contador = 1;
  
  A_artificial = [A eye(m)];
  
  indB_original = [1:n]';
  indNB = [1:n]';
  
  indB = setdiff(1:(m+n), indNB)';
  indB_estatico = indB;
  
  Binv = eye(m);
  
  c_artificial = [zeros(1, n) ones(1, m)]';
  
# Na primeira iteração x(indB) = b para as variáveis básicas y_1 ... y_m
  x = zeros(1, size(c_artificial, 1))';
  x(indB) = b;
   
  #inicializa uma matriz para guardar os valores do vetor u. Para aplicar o algoritmo do Bertsimas
  matriz_u = zeros(m, n);
 
	# Loop principal do programa
	while contador == 1  || min(c_barra) < 0
		
   
		# Primeira iteração
		fprintf('\n')	
		fprintf('***************** ITERANDO %d *****************\n', contador-1)
		fprintf('\n')	
		
    # variáveis não básicas
    var_nb = setdiff(1:length(x), indB)';
    
    funcao_custo = - sum(x(intersect(indB_estatico, indB)));
    

    fprintf('Índices e valores das variáveis básicas:\n')

		for i = 1: size(indB,1)
			fprintf('%d   %f\n',indB(i),x(indB)(i))
		endfor
		fprintf('\n')
    
    # Definindo os vetores associados às variáveis básicas
		# Matriz básica B 
		B = A_artificial(:,indB);
    
    fprintf('Matriz inversa da base:\n')
    disp(Binv)
    fprintf('\n')
    
		# Variáveis básicas associadas à solução básica inicial
		x_b = x(indB);

		# Vetor de custo associado às variáveis básicas
		c_b = c_artificial(indB);

		# Cálculo do vetor p = (c_B)^T * B^-1
    
    
		fprintf('Vetor p:\n')

		p = c_b' * Binv;
    
    
		for i = 1: size(p,2)
			fprintf("   %f\n", p(i))  
		endfor
		fprintf('\n')

		# Inicializa o vetor c_barra com zeros
		c_barra = zeros(1, n)';
    
		# Faz o cálculo de vetor de custos c_j_barra = c_j - (p * A_j); em que p = (c_B)^T * B^-1
		for j = var_nb'
      # Cálculo do custo reduzido para as variáveis não básicas
      c_barra(j) = c_artificial(j) - p * A_artificial(:,j);
		endfor
    
    
		# Printa os custos reduzidos calculados para as variáveis não básicas, uma vez que para as variáveis básicas este custo é 0 
		fprintf('Custos Reduzidos:\n')

		# Vetor com as variáveis não básicas
    
    # faz a diferença de conjuntos para determinar os índices das variáveis não básicas
		for i = var_nb'
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
		u(indB) = Binv * A_artificial(:,j);
    
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
    indice_theta_minimo = find(theta == theta_minimo)(1);
    
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
     
    # Índices básicos: j enrra na base e indice_theta_minimo sai da base (l)
		

	  # Faz o cálculo da nova solução viável básica para as variáveis básicas
    
    # Atribui theta_estrela à posição j de x
    
    # Índices básicos: j entra na base e indice_theta_minimo sai da base (l)
		
    x(j) = theta_minimo;
   
    # Faz o cálculo da nova solução viável básica para as variáveis básicas
		x(indB) = x_b - theta_minimo*u(indB);
 
		# T é a matriz aumentada [inv(B) , vetor u]
		T = [Binv u(indB)];
    
    
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
    
    # salva o vetor u
    matriz_u(1:end, j) = T(:, end);
	
		# Aumenta o contador
		contador = contador + 1;
		
	endwhile
	
# se tiver solução (ind = 0), retorna o custo ótimo, os índices e valores das variáveis básicas	
  if ind == 0
        
    if funcao_custo == 0
      # verificar se tem variaveis auxiliares no tableau (variavel indB), se tiver tirar ela seguindo (5)
      
      # faz a intersecao para verificar se existem variáveis auxiliares nas variáveis básicas atuais
        intersecao = intersect(indB, indB_estatico);
        
      # se a matriz u não estiver completa, completa com as vars. que faltam:
        for j = setdiff(indB_original, indB)
           matriz_u(1:end, j) = Binv * A(:, j);  
        endfor
       
        # variável matriz u estendida com Binv
        matriz_u = [matriz_u Binv];
     
      
      # para cada uma dessas variáveis, performa (5) - página 117 Bertsimas
        for i = intersecao
        # a linha i é redundante
          
          linha = find(indB_estatico == i);   
        
          #matriz_original = 
        # se a soma dos módulos for 0, então não existem elementos != 0 e a variável l (x_l) é eliminada do problema
          if sum(abs(matriz_u(linha, indNB))) == 0
            
                    
            fprintf('Variável redundante. Elimina-se a variável: %d\n', indB(linha)) 
            fprintf('\n')
            
            # elimina a linha na matriz e o índice no vetor dos índices básicos
            matriz_u(linha, :) = [];
            indB(linha) = [];
            
         
          # se a soma não for 0, então pelo menos um elemento é diferente de 0, para ele, fazemos:
          else
      
            # para cada variável que não deve estar na base original (l-ésima linha), faz o processo de tirar da base
            
           
              
             indice_nao_zero = find(matriz_u(linha, indB_original) != 0)(1);
            
             fprintf('Entra na base: %d\n', indice_nao_zero)
             fprintf('\n')
    
    
             fprintf('Sai da base: %d\n', indB(linha)) 
             fprintf('\n')
     
         
             matriz_u(linha, :) = matriz_u(linha, :)/matriz_u(linha, indice_nao_zero);
          
          
         
             for j = 1:m
               if j != indice_nao_zero
                 matriz_u(j, :) = - matriz_u(j, indice_nao_zero) * matriz_u(linha, :) + matriz_u(j, :);
               endif
             endfor   
          
             
             indB(linha) = indice_nao_zero;
                
          endif
        
        endfor
        
        fprintf('Solução viável básica encontrada:\n')
	      fprintf('\n')
	      fprintf('Índices e solução ótima:\n')
	      
        for i = 1:size(indB_original, 1)
		      fprintf('%d   %f\n', i, x(i))  
	      endfor
	      
    
        A = matriz_u(:, 1:n);
        x = x(1:n);
        c_barra = zeros(n,1);
        
        m = size(A, 1);
        n = size(A, 2);
        

      
        [ind v] = simplex_revisado(A,b,c,m,n,x,indB, eye(size(indB, 1)));
        
        # se o problema for ilimitado, retorna x vazio, vetor direção e variável ind
        if ind == -1
          ind;
          x = [];
          v = v;
        endif
    endif
      
    # Verificando a função custo é > 0 nas variáveis aux., se for, o problema é inviável e o algoritmo acaba
    if funcao_custo != 0
        ind = 1;
        v = [];
        x = [];
        fprintf('Problema inviável. Fim de programa.\n')
	      fprintf('\n')
    endif
  
  endif
  
 endfunction
