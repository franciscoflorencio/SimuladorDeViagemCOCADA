---
# IMPORTANTE!!! FAZER DOWNLOAD DO DATASET: https://drive.google.com/file/d/19g0f6WNUiDdSOyA1Ef_tK8Xcl2tjI9aN/view?usp=sharing
---
# Simulador de Viagem entre Capitais Brasileiras

## Introdução

Imagine que você é uma pessoa planejando viajar de ônibus entre capitais brasileiras. Para otimizar sua viagem, você deseja encontrar o trajeto mais rápido entre duas cidades, considerando apenas rotas possíveis entre estados que compartilham fronteiras. Este projeto simula exatamente essa situação!

A ideia surgiu da necessidade de planejar viagens mais eficientes e seguras, considerando os dados reais de distâncias rodoviárias brasileiras disponíveis no site [Brazil Road Distances](https://rfsaldanha.github.io/data-projects/brazil_road_distances.html) de Rodrigo Saldanha. Utilizando esses dados, foi criado um algoritmo capaz de determinar o menor caminho entre capitais, utilizando conceitos de álgebra linear e estruturas algébricas em Julia.

## O Contexto da Viagem

João, nosso protagonista, quer viajar de Rio Branco (AC) até Maceió (AL) de ônibus. Porém, ele percebe que nem todas as capitais estão diretamente conectadas. João precisa planejar cuidadosamente para garantir que sua rota seja viável e que ela chegue ao destino no menor tempo possível. Para isso, utilizamos um modelo matemático que leva em conta:
- As distâncias reais entre as capitais.
- As conexões entre estados vizinhos.
- A rota com a menor distância acumulada.

O algoritmo descrito neste relatório permite que João insira sua cidade de origem e destino e, em poucos segundos, obtenha a melhor rota possível.

## Conceitos Fundamentais

### Monóides e Semianéis

Um monóide é definido como uma 3-tupla $(S, \oplus, 0)$, onde:
- $S$ é um conjunto não-vazio;
- $\oplus$ é uma operação binária chamada "soma";
- $0$ é o elemento neutro em relação a $\oplus$, tal que $x \oplus 0 = 0 \oplus x = x$.

O monóide satisfaz as propriedades de totalidade, associatividade e existência de elemento neutro.

Um semianel é uma generalização do monóide, definido como uma 5-tupla $(S, \oplus, \cdot, 0, 1)$, com:
- $(S, \oplus, 0)$ sendo um monóide comutativo;
- $(S, \cdot, 1)$ formando outro monóide;
- A operação $\cdot$ distribuindo sobre $\oplus$.

### Fecho e Ponto Fixo

Um semianel fechado é uma extensão de um semianel $(S, \oplus, \cdot, 0, 1)$ com uma operação adicional de fecho, denotada por $*$. Formalmente, o fecho de $x \in S$, denotado por $x^*$, é definido como:
\[ x^* = 1 \oplus x \oplus (x \cdot x) \oplus (x \cdot x \cdot x) \oplus \cdots \]

A operação de fecho é útil em várias aplicações, como resolver sistemas lineares em álgebra tropical e calcular caminhos mínimos em grafos.

### Aplicação em Grafos

No caso de grafos ponderados, o fecho é utilizado para encontrar o menor caminho entre os vértices. Especificamente, no semianel min-plus, o fecho de uma matriz de adjacência $G$ é definido como:
\[ G^* = \bigoplus_{n=1}^\infty G^n \]

A convergência para o ponto fixo ocorre quando o valor de $G^{n+1}$ não difere de $G^n$.

## Explicação Detalhada do Código

### Configuração Inicial

O código inicia com as importações necessárias:
```julia
using DataFrames
using CSV
using LinearAlgebra
```

### Estrutura `Path`

A estrutura `Path` é utilizada para armazenar trajetos entre cidades:
```julia
struct Path
    p::Union{Nothing, Vector{Int}}
    d::Float64
end
```

### Métodos Personalizados para `Path`

Redefinição de operadores para comparar e concatenar caminhos:
```julia
function Base.:+(a::Path, b::Path)
    isnothing(a.p) && return b
    isnothing(b.p) && return a
    return a.d < b.d ? a : b
end

function Base.:*(a::Path, b::Path)
    isnothing(a.p) && return Path(nothing, Inf)
    isnothing(b.p) && return Path(nothing, Inf)
    @assert last(a.p) == first(b.p)
    return Path([a.p; b.p[2:end]], a.d + b.d)
end

Base.zero(::Type{Path}) = Path(nothing, Inf)
```

### Matriz Identidade

Criação de uma matriz identidade para `Path`:
```julia
function Identity(dim, tipo)
    return [Path(i == j ? [i] : nothing, i == j ? 0.0 : Inf) for i in 1:dim, j in 1:dim]
end
```

### Multiplicação de Matrizes

Realização da multiplicação de duas matrizes de caminhos:
```julia
function Base.:*(A::Matrix{Path}, B::Matrix{Path})
    m, n = size(A)
    p = size(B, 2)
    result = Matrix{Path}(undef, m, p)
    
    for i in 1:m
        for j in 1:p
            paths = [A[i,k] * B[k,j] for k in 1:n]
            result[i,j] = reduce(+, filter(path -> path.p !== nothing, paths), init=Path(nothing, Inf))
        end
    end
    
    return result
end
```

### Construção do Grafo de Capitais

Carregamento e filtragem dos dados:
```julia
df = CSV.read("dist_brasil_clean.csv", DataFrame)
df_fronteiras = filter(row -> (row.orig, row.dest) in fronteiras || (row.dest, row.orig) in fronteiras, df)
```

Construção da matriz de adjacência:
```julia
matriz_adj = Identity(n, Path)

for row in eachrow(df_fronteiras)
    i = capitais_dict[row.orig]
    j = capitais_dict[row.dest]
    matriz_adj[i, j] = Path([i, j], row.dist)
    matriz_adj[j, i] = Path([j, i], row.dist)
end
```

### Algoritmo de Menor Caminho Utilizando o Método do Ponto Fixo

Implementação do algoritmo de menor caminho:
```julia
function menor_caminho(matriz_adj, origem, destino)
    n = size(matriz_adj, 1)
    pm = copy(matriz_adj)
    fixed_points = 0

    while fixed_points < 2
        pm_old = copy(pm)
        pm = pm * pm
        if pm == pm_old
            fixed_points += 1
        else
            fixed_points = 0
        end
    end

    return pm[origem, destino]
end
```

### Consulta Interativa

Consulta interativa para obter o menor caminho:
```julia
println("Digite o índice da cidade de partida:")
origem_idx = parse(Int, readline())
println("Digite o índice da cidade de destino:")
destino_idx = parse(Int, readline())

caminho = menor_caminho(matriz_adj, origem_idx, destino_idx)
println("O menor caminho entre $origem_idx e $destino_idx é $caminho")
```

## Mapeamento de Capitais e Fronteiras

### Representação Geográfica com Códigos Municipais do IBGE

Exemplo de códigos municipais das capitais brasileiras:
- Rio Branco: 1200401
- Maceió: 2704302
- Belo Horizonte: 3106200

Construção das fronteiras no algoritmo:
```julia
fronteiras = [
    (1100205, 5103403),  # Porto Velho (RO) e Cuiabá (MT)
    (1200401, 1302603),  # Rio Branco (AC) e Manaus (AM)
    (3106200, 3304557),  # Belo Horizonte (MG) e Rio de Janeiro (RJ)
    # Outras fronteiras definidas manualmente...
]
```

## Conclusão

Esse trabalho trouxe um resultado satisfatório e interessante, algo bem legal de ver acontecendo e se observar. Agradecimentos especiais ao professor João Paixão e ao monitor de disciplina Matheus do Ó pelo suporte e orientação ao longo do projeto.

## Referências

- [Distances in Julia](https://bkamins.github.io/julialang/2022/10/21/distances.html)
- [Brazil Road Distances](https://rfsaldanha.github.io/data-projects/brazil_road_distances.html)
- [Algoritmos Genéricos para Multiplicação Matricial](https://github.com/Mth0/Algoritmos-Genericos-Multiplicacao-Matricial)
- Fun with semirings: A functional pearl on the abuse of linear algebra, S. Dolan, Proceedings of the 18th ACM SIGPLAN International Conference on Functional Programming, 2013.
- CFG Parsing and Boolean Matrix Multiplication, F. C. Ebert, 2007.
- Parsing Techniques - A Practical Guide, D. Grune, C. J. H. Jacobs, Springer, 2008.
- What is the shortest path to AlphaTensor?, B. Kamiński.
- Algebraic Structures for Transitive Closure, D. J. Lehmann, Theoretical Computer Science, 1977.

---
