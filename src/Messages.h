/*****************************************************************************/
/***																	   ***/
/***		Modificado por G.J.L.P.										   ***/
/***		Añadidos nuevos mensajes que indican falta de algún 	 	   ***/
/***		Fichero de Configuración (No específico para ningún	  		   ***/
/***		problema) o nuevos errores.				  					   ***/
/***																	   ***/
/*****************************************************************************/

#ifndef RLFAP_MESSAGES
#define RLFAP_MESSAGES

#ifndef MAX_BUFFER
#define MAX_BUFFER 200
#endif

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

inline void show_message(int value)
{
	switch (value)
	{
		case 1: std::cout << std::endl << "Error: number of arguments in the execution call is incorrect !!"
			     << std::endl; break;
		case 2: std::cout << std::endl << "Error: It's imposible find Configuration file !!" << std::endl;
			break;
		/* Específicos de RLFAP */
		case 3:	std::cout << std::endl << "Error: It is imposible find the Celar problem definition file (cst.txt) !!"
			     << std::endl; break;
		case 4: std::cout << std::endl << "Error: It is imposible find the Celar domains file (dom.txt) !!"
			     << std::endl; break;
		case 5:	std::cout << std::endl << "Error: It is imposible find the Celar links file (var.txt) !!"
			     << std::endl; break;
		case 6:	std::cout << std::endl << "Error: It is imposible find the Celar constraints file (ctr.txt) !!"
			     << std::endl; break;
		/* Fallos de Memoria */
		case 7:	std::cout << std::endl << "Error: No avalible memory for \"malloc\" operation !!"  << std::endl;
			break;
		case 8:	std::cout << std::endl << "Error: in \"free\" operation !!"  << std::endl;
			break;
		/* Específicos del MaxCut */
		case 9:	std::cout << std::endl << "Error: It is imposible find the Maxcut file (Maxcut.txt) !!"
			     << std::endl; break;
		/* Genéricos de Falta de ficheros de configuracion  adicionales al mensaje 2 */
		case 10: std::cout << std::endl << "Error: It's imposible find Configuration file (Config.cfg) !!"
			      << std::endl; break;
		case 11: std::cout << std::endl << "Error: It's imposible find Skeleton Configuration File (Ske.cfg) !!"
			      << std::endl; break;
		case 12: std::cout << std::endl << "Error: It's imposible find Instance Problem File !!" << std::endl;
			 break;
		case 13: std::cout << std::endl << "Error: It's imposible find Resultate File !!" << std::endl;
			 break;
		case 14: std::cout << std::endl << "Error: Index out of Range !!" << std::endl;
			 break;
		default: std::cout << std::endl << "Unkown Error !!" << std::endl;
	}

	std::cout << std::endl << " " << std::endl;
	exit(-1);
}

inline void continue_question()
{
	fflush(stdout);
	std::cout << std::endl << "Press any key to continue..." << std::endl;
	fflush(stdin);
	getc(stdin);
}

inline void get_path(const char *source,char *target)
{
	int last = 0;

	for(int i = 0; i < strlen(source); i++)
	{
		target[i] = source[i];
		if(target[i] == '/')
			last = i;
	}
	target[last+1] = '\0';
}

inline unsigned count_lines(char *file_name) // returns the number of lines of a file
{
	char line[MAX_BUFFER];
	FILE *file;
	int count=0;

	if ((file=fopen(file_name,"r"))==NULL)
	{
		fflush(stdout);
		printf("File not found !");
	}
  
	while (!feof(file))
	{
		if (fgets(line,MAX_BUFFER,file)) count++;    		     	
		else
		{
			fclose(file);
			break;
		}
	}	
	return count;
}

#endif
