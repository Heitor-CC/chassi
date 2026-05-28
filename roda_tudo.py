# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
import job
import sys
import os

# muda a pasta de trabalho do abaqus
caminho_pasta = r"C:\Users\Heitor\Documents\Github\chassi"
if caminho_pasta not in sys.path:
    sys.path.append(caminho_pasta)

# importa o modelador de chassi =====> se a linha "from constroi_chassi_oficial import GerarModeloChassi, modelar" não estiver aqui, coloque!!!
from constroi_chassi_oficial import GerarModeloChassi, modelar


# diametros comuns de tubos de pvc. O utilizado oficialmente será o de 25mm
diametros = [3.2, 2.5, 2.0]

dados = {
    3.2: {'comp_joelho': 2.10, 'cota_z': 1.80},   # bitola de 32mm
    2.5: {'comp_joelho': 1.69, 'cota_z': 1.44},   # bitola de 25mm
    2.0: {'comp_joelho': 1.40, 'cota_z': 1.14}    # bitola de 20mm
}


# medidas do chassi padrão
L_fixo = 53.5
h_fixo = 15.0
b_fixo = 41.0
t_fixo = 0.15
R_roda = 11.0
peso_crianca = 50.0
prop_banco = 0.5


for d in diametros:
    print('--> Gerando modelo para o diametro: ' + str(d*10) + ' mm')

    nome_job = 'Simulacao_diametro_' + str(int(d*10)) + 'mm'
    nome_modelo = 'Chassi-'+str(int(d*10))+'mm'

    comp_joelho = dados[d]['comp_joelho']
    cota_z = dados[d]['cota_z']

    # cria o chassi com o novo diametro
    chassi_atual = GerarModeloChassi(
        L=L_fixo, h=h_fixo, b=b_fixo, D=d, t=t_fixo,
        R=R_roda, peso=peso_crianca, proporcao_banco=prop_banco, comp_joelho_cm=comp_joelho, nome=nome_modelo
    )

    # modela
    modelar(chassi_atual)

    # imprime as barras para cada diametro
    print('lista de corte para bancada:')
    chassi_atual.gerar_lista_de_corte(cota_z)

    # criação do job
    mdb.Job(name=nome_job, model=nome_modelo,
            description='Chassi com diametro ' + str(int(d*10)) + 'mm')

    # submissão
    print('    Submetendo Job: ' + nome_job)
    mdb.jobs[nome_job].submit()


print('------------------------------------------------')
print('Sucesso! Todas as analises de diametro foram concluidas.')
