function [ind x v] = simplex_tab(tableau)
  
  fprintf('\n*********************** Simplex: Fase 1 -  Método Tableau ************************\n\n')  
  
  x = tableau(2:end, 1);
  
  # faz a regra para ter b >= 0
  ind_negativo_b = find(x < 0);
  
  # multiplica as linhas do tableau por -1 nos índices ind_negativo_b
  tableau(ind_negativo_b + 1, :) = - tableau(ind_negativo_b + 1, :); 
  
  A = tableau(2:end, 2:end);
  m = size(tableau(2:end, 1), 1);
  n = size(tableau(1:end, 2:end), 2);
  x = tableau(2:end, 1);
  
  custo_original = tableau(1, 2:end);
  
  indNB = [1:n];
  indB = setdiff(1:(m+n), indNB);
  indB_estatico = indB;
  
  
  A_artificial = [A eye(m)];

  # Passando para o formato artificial
  custo_reduzido = - [ones(1, m) * A_artificial(:, indNB) zeros(1, m)];
  funcao_custo_artificial = - sum(x);
  
  tableau = [funcao_custo_artificial custo_reduzido; [x A_artificial]];
    
  contador = 1;
  
  # while loop para a realização dos cálculos
  while contador == 1  || min(c_barra) < 0
    
    funcao_custo = tableau(1, 1); # valor função custo
    
    x = tableau(2:end, 1); # vetor x
    
    c_barra = tableau(1, 2:end); # custos reduzidos
    
  
  # Primeira iteração
    fprintf('\n')	
    fprintf('********************************* ITERANDO %d *********************************\n', contador-1)
    fprintf('\n')	
    
    # A entrada do tableau deverá ser feita por linhas, ou seja, uma linha inteira de cada vez,
    # assim este programa funcionará corretamente sem problemas de formato.
    
    
    fprintf('Tableau:\n')
    disp(tableau)
    fprintf('\n')

    # Verificando o vetor de custos reduzidos
		
    if min(round(c_barra)) >= 0
      ind = 0;
			break  
    endif
    
    # Performando o Simplex
    # aplica a regra de Bland (regra do menor índice)
    # se tiver algum índice negativo, ele será o índice j da direção básica viável
    j = find(c_barra < 0 , 1);
    
    # Construção de u (1 x n), o vetor direção (j-ésima direção básica)
    u = tableau(2:end,j+1);
    
    
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
    
    
    # se a funcao_custo for < 0, performa o Simplex até se obter a funcao_custo > 0 ou = 0

    # Uma solução viável foi encontrada para o problema original
    
    # aumenta o contador
    contador = contador + 1;
  
  endwhile
 
  
  if ind == 0
    
    funcao_custo = round(funcao_custo);

    if funcao_custo == 0
    # verificar se tem variaveis auxiliares no tableau (variavel indB), se tiver tirar ela seguindo (5)
      
    # faz a intersecao para verificar se existem variáveis auxiliares nas variáveis básicas atuais
      intersecao = intersect(indB, indB_estatico);
    
      # para cada uma dessas variáveis, performa (5) - página 117 Bertsimas
      for i = intersecao
        # a linha i é redundante
        
        linha = find(indB_estatico == i);
        
        # se a soma dos módulos for 0, então não existem elementos != 0 e a variável i (x_i) é eliminada do problema
        if sum(abs(tableau(linha+1, indNB+1))) == 0
          
          fprintf('Variável redundante. Elimina-se: %d\n', indB(linha)) 
          fprintf('\n')
         
          tableau(linha+1, :) = [];
          indB(linha) = [];
          x = tableau(2:end, 1);
         
          # se a soma não for 0, então pelo menos um elemento é diferente de 0, para ele, fazemos:
        else
          
          indice_nao_zero = find(tableau(linha+1, indNB+1) != 0)(1);
            
          fprintf('Entra na base: %d\n', indice_nao_zero)
          fprintf('\n')
            
          fprintf('Sai da base: %d\n', indB(linha)) 
          fprintf('\n')

          tableau(linha+1, :) = tableau(linha+1, :)/tableau(linha+1, indice_nao_zero+1);
           
          for j = 1:m+1
            if j != indice_nao_zero+1 && tableau(j, indice_nao_zero+1) != 0
              tableau(j, :) = - tableau(j, indice_nao_zero+1) * tableau(linha+1, :) + tableau(j, :);
            endif
          endfor  
          
         
         
          indB(linha) = indice_nao_zero;
        endif     
      endfor
      
      fprintf('Solução viável básica encontrada:\n')
	    fprintf('\n')
	    fprintf('Índices e solução ótima:\n')
      
      x_final = zeros(n, 1);
      x_final(indB) = x;
    
      for i = 1:size(x_final, 1)
       fprintf("%d   %f\n", i, x_final(i))
      endfor

      # chamando a Fase 2 do Método Simplex Tableau
      for i = setdiff(1:n, indB)
        tableau(1, i+1) = custo_original(i) - custo_original(indB) * tableau(2:end, i+1);
      endfor
      
      # calcula o valor da função custo, por definição: - c_{B}^{T} * B^{-1} b = - c_{B}^{T} * x
      tableau(1,1) = - dot(custo_original(indB), x_final(indB));
      
      x = x_final;
      tableau = tableau(1:end, 1:(n+1));
      
      [ind v] = simplex_tableau(indB, tableau);
      
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