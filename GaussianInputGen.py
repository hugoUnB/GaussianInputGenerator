import os
import re
class ChalcInputGen:
    periodic_table = {
        "1": "H", "2": "He",
        "3": "Li", "4": "Be", "5": "B", "6": "C", "7": "N", "8": "O", "9": "F", "10": "Ne",
        "11": "Na", "12": "Mg", "13": "Al", "14": "Si", "15": "P", "16": "S", "17": "Cl", "18": "Ar",
        "19": "K", "20": "Ca", "21": "Sc", "22": "Ti", "23": "V", "24": "Cr", "25": "Mn", "26": "Fe", "27": "Co",
        "28": "Ni", "29": "Cu", "30": "Zn", "31": "Ga", "32": "Ge", "33": "As", "34": "Se", "35": "Br", "36": "Kr",
        "37": "Rb", "38": "Sr", "39": "Y", "40": "Zr", "41": "Nb", "42": "Mo", "43": "Tc", "44": "Ru", "45": "Rh",
        "46": "Pd", "47": "Ag", "48": "Cd", "49": "In", "50": "Sn", "51": "Sb", "52": "Te", "53": "I", "54": "Xe"
    }

    def __init__(self, path=os.getcwd()):
        """
        Esta classe extrai as coordenadas moleculares de uma série de arquivos 'gjf'
        e cria inputs personalizados a partir destas informações

        :param path: 'String' do diretório que contém os arquivos 'gjf'
        """
        self.pathin = path
        files_gjf = [file for file in os.listdir(path) if '.gjf' in file]
        files_out = [file for file in os.listdir(path) if '.out' in file or '.log' in file]

        self.coord = {}
        for file in files_gjf:
            with open(os.path.join(path,file), 'r') as txt:
                self.coord[file] = self.__extract_coord_gjf(txt.readlines())

        for file in files_out:
            with open(os.path.join(path,file), 'r') as txt:
                self.coord[file] = self.__extract_coord_out(txt.readlines())

    def __extract_coord_gjf(self,txt):
        """
        Este método extrai as linhas que contém as coordenadas de um arquivo 'gjf'
        :param txt: 'Lista' contendo as linhas do arquivo gjf
        :return:  'Lista' contendo apenas as linhas referentes as coordenadas moleculares
        """
        begin = 0
        for i, line in enumerate(txt):
            result = re.match(r'^[+-]?[0-9]?\s?[0-9]$', line)
            if result != None:
                begin = i
                break
        if begin == 0:
            return None
        for i, line in enumerate(txt[begin:]):
            result = re.match(r'^\n$', line)
            if result != None:
                end = i + begin
                break
        return txt[begin:end]

    def __extract_coord_out(self,txt):
        """
        Este método extrai as coordenadas de um arquivo 'out' ou 'log' já no formato de escrita do 'gjf'
        :param txt: 'Lista' contendo as linhas do arquivo 'out' ou 'log'
        :return:  'Lista' contendo apenas as linhas referentes as coordenadas moleculares
        """
        idx = 0
        for i, line in enumerate(txt):
            result = re.search('Standard orientation', line)
            if result != None:
                idx = i + 5
        if idx == 0: return None

        coord = []
        for line in txt[idx:]:
            if re.search('----------------------------', line) != None:
                break
            ln = line.split()
            del ln[0]
            del ln[1]
            ln[0] = self.periodic_table[ln[0]]
            ln.append('\n')
            ln = '          '.join(ln)
            coord.append(ln)
        return coord

    def __gjf_head(self,name,storage,proc,mem,calc_line):
        """
        Este método cria o cabeçalho de um arquivo gjf
        :param name: 'String' que será o título da molecula
        :param storage: 'String' diretório que o usuário quer utilizar dentro de /Storage01/
        :param proc:  'Float' qtd de processadores requeridos
        :param mem: 'Float' qtd de memória requerida
        :param calc_line: 'String' linha de comandos para o cálculo no Gaussian
        :return: 'Lista' contendo as linhas de um cabeçalho de um arquivo gjf
        """
        return [
            '%rwf=/Storage01/{}\n'.format(storage),
            '%int=/Storage01/{}\n'.format(storage),
            '%d2e=/Storage01/{}\n'.format(storage),
            '%nosave\n',
            '%nprocshared={}\n'.format(proc),
            '%mem={}GB\n'.format(mem),
            '%chk={}\n'.format(name),
            '# {} \n'.format(calc_line),
            '\n',
            '{}\n'.format(name),
            '\n']

    def write_gjf(self,
                  path = 'inputs',
                  deut = True,
                  storage = 'hugo/chalcona',
                  proc = 8,
                  mem = 16,
                  calc_line = 'opt freq=noraman X3LYP/6-31+G(d) volume scrf=(smd,solvent=water) scf=xqc',
                  name_add = ''):
        """
        Este método gera os arquivos gjf personalizados
        :param path: 'String' diretório onde serão criados os gjf
        :param deut: 'Boolean' use False se quiser transformar todos os deutérios em hidrogenios
        :param storage: 'String' diretório que o usuário quer utilizar dentro de /Storage01/
        :param proc: 'Float' qtd de processadores requeridos
        :param mem: 'Float' qtd de memória requerida
        :param calc_line: 'String' linha de comandos para o cálculo no Gaussian
        """
        self.pathout = os.path.join(self.pathin, path)
        try:os.mkdir(self.pathout)
        except:pass

        if deut == False:
            for file in self.coord:
                for i, line in enumerate(self.coord[file]):
                    if '(Iso=2)' in line:
                        self.coord[file][i] = line.replace('(Iso=2)','       ')
        self.name_add = name_add
        for file in self.coord:
            idx = file.index('.')
            filename = self.__change_name(file)

            with open(os.path.join(self.pathout,filename),'w') as w:
                w.writelines(self.__gjf_head(file[:idx],storage,proc,mem,calc_line))
                w.writelines(self.coord[file])
                w.write('\n\n')

    def __pbs_head(self,name,nodes,mem,walltime,g0):
        """
        Este método cria o cabeçalho personalizado de um arquivo pbs
        :param name: 'String' nome que aparecerá na fila do cluster
        :param nodes: 'Float' qtd de processadores requeridos
        :param mem: 'Float' qtd de memória requerida
        :param walltime: 'String' ou 'Float' número de horas de cálculo para entrar na fila
        :param g0: 'Float' versão do gaussian para carregamento do módulo correto
        :return: 'Lista' contendo as linhas do cabeçalho de um arquivo arquivo pbs
        """
        return [
            '##\n',
            '##\n',
            '#PBS -S /bin/bash\n',
            '#PBS -l nodes=1:ppn={}\n'.format(nodes),
            '#PBS -l mem={}GB\n'.format(mem),
            '#PBS -l walltime={}:00:00\n'.format(walltime),
            '#PBS -N {}\n\n'.format(name),
            'cd $PBS_O_WORKDIR\n\n',
            'echo "inicio do job: "`date`\n',
            'echo "Hostname: " `hostname`\n',
            'echo "PWD: "$PWD\n\n',
            'module load softwares/gaussian-{}\n\n'.format({9:r'09/pgi', 16:r'16/b01'}[g0])]

    def __change_name(self,name):
        """
        Este método altera o nome do arquivo
        :param name: 'String' nome de entrada
        :return: 'String' nome de saída
        """
        idx = name.index('.')
        if '.gjf' not in name:
            filename = '{}_{}.gjf'.format(name[:idx], self.name_add)
        elif self.name_add != '':
            filename = '{}_{}.gjf'.format(name[:idx], self.name_add)
        else:
            filename = name
        return filename

    def write_pbs(self,
                  name = 'HugoChalc',
                  nodes = 8,
                  mem = 16,
                  walltime = 1440,
                  g0 = 16):
        """
        Este método cria o arquivo PBS para submeter os inputs criados no cluster
        :param name: 'String' nome que aparecerá na fila do cluster
        :param nodes: 'Float' qtd de processadores requeridos
        :param mem: 'Float' qtd de memória requerida
        :param walltime: 'String' ou 'Float' número de horas de cálculo para entrar na fila
        :param g0: 'Float' versão do gaussian para carregamento do módulo correto (suporta apenas 9 e 16)
        """
        with open(os.path.join(self.pathout,'fila.pbs'),'w') as w:
            w.writelines(self.__pbs_head(name,nodes,mem,walltime,g0))

            for file in self.coord:
                idx = file.index('.')
                filename = self.__change_name(file)
                w.write('g{} {}\n'.format({9:'09', 16:'16'}[g0],filename))
                w.write('sleep 3\n\n')


if __name__ == '__main__':

    a = ChalcInputGen('example_files')
    a.write_gjf(path='saida_normal', deut=False, name_add='normal')
    a.write_pbs()

    b = ChalcInputGen('example_files')
    b.write_gjf(path='saida_deut', deut=True, name_add='deut', storage='hugo/chalcona_deut')
    b.write_pbs('HugoChalc_deut')

