/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: misc.cc,v 1.22 2003-06-11 21:33:18 krueger Exp $ */
/*
* ABSTRACT: lib parsing
*/

#include <stdlib.h>
#include <mod2.h>
#include <febase.h>
#include <grammar.h>
#include <ipid.h>
#include <ipshell.h>
#include <structs.h>
#include <subexpr.h>
#include <tok.h>
#include <regex.h>

#include "modgen.h"
#include "typmap.h"

#define NEW_PARAM 1

#define SYSTYP_NONE   0
#define SYSTYP_LINUX  1
#define SYSTYP_HPUX9  2
#define SYSTYP_HPUX10 3

#if 1
#  define logx printf
#else
#  define logx
#endif

char *DYNAinclude[] = {
   "",
   "#include <dlfcn.h>",
   "#include <dl.h>",
   "#include <dl.h>",
   "#include <>",
   "#include <>",
   "#include <>",
   NULL
};

int systyp = SYSTYP_NONE;

BOOLEAN    expected_parms;
int        cmdtok;
BOOLEAN siq=FALSE;
char *lastreserved=NULL;
#define SELF_CMD MAX_TOK+1

extern void enter_id(FILE *fp, idtyp t, char *name, char *value,
                     int lineno, char *file);
extern void write_enter_id(FILE *fp);
extern void write_add_singular_proc(FILE *fp);

static void  mod_write_ctext(FILE *fp, FILE *fp_in);
char type_conv[MAX_TOK][32];

extern void write_intro(moddefv module);

/*=============== types =====================*/
struct sValCmdTab
{
  short cmd;
  short start;
};

#include "../../Singular/iparith.inc"

/*=================== general utilities ============================*/
int IsCmd(char *n, int & tok)
{
  int an=1;
  int i,v;
  int en=LAST_IDENTIFIER;
  
  if( strcmp(n, "SELF") == 0) {
    tok = SELF_CMD;
    logx("IsCmd: [%d] %s\n", tok, n);
    return tok;
  }
  
  if( strcmp(n, "none") == 0) {
    tok = NONE;
    logx("IsCmd: [%d] %s\n", tok, n);
    return tok;
  }
  
  loop
  {
    if(an>=en-1)
    {
      if (strcmp(n, cmds[an].name) == 0)
      {
        i=an;
        break;
      }
      else if (strcmp(n, cmds[en].name) == 0)
      {
        i=en;
        break;
      }
      else
      {
        return 0;
      }
    }
    i=(an+en)/2;
    v=strcmp(n,cmds[i].name);
    if(v<0)
    {
      en=i;
    }
    else if(v>0)
    {
      an=i;
    }
    else /*v==0*/
    {
      break;
    }
  }
  lastreserved=cmds[i].name;
  tok=cmds[i].tokval;
  if(cmds[i].alias==2)
  {
    if(trace)printf("outdated identifier `%s` used - please change your code",
    cmds[i].name);
    cmds[i].alias=1;
  }
#if 0
  if (!expected_parms)
  {
    switch (tok)
    {
      case IDEAL_CMD:
      case INT_CMD:
      case INTVEC_CMD:
      case MAP_CMD:
      case MATRIX_CMD:
      case MODUL_CMD:
      case POLY_CMD:
      case PROC_CMD:
      case RING_CMD:
      case STRING_CMD:
        cmdtok = tok;
        break;
    }
  }
#endif
  logx("IsCmd: [%d] %s\n", tok, n);
  
  if( (cmds[i].toktype==ROOT_DECL) ||
      (cmds[i].toktype==ROOT_DECL_LIST) ||
      (cmds[i].toktype==RING_DECL) ||
      (cmds[i].toktype==IDEAL_CMD) ||
      (cmds[i].toktype==INTMAT_CMD) ||
      (cmds[i].toktype==MODUL_CMD) ||
      (cmds[i].toktype==MATRIX_CMD))// ||
//      ((cmds[i].toktype>=DRING_CMD) && (cmds[i].toktype<=VECTOR_CMD))) 
    return cmds[i].toktype;
  return 0;
}

char * decl2str(int n, char *name)
{
  switch(n)
    {
#   include "decl.inc"

    /* first and last entry of tok.h cannot be grepped */
      //case ATTRIB_CMD: strcpy(name,"ATTRIB_CMD"); break;
    case POLY_CMD: strcpy(name,"POLY_CMD");       break;
    case MATRIX_CMD: strcpy(name,"MATRIX_CMD");       break;
    case MAX_TOK: strcpy(name,"MAX_TOK");       break;
    default: strcpy(name,"(null)");
  }
  return(name);
}

/*========================================================================*/
struct valid_cmds_def 
{
  char *name;
  void (*write_cmd)(moddefv module, procdefv pi, void *arg = NULL);
  cmd_token id;
  cmd_type  type;
  int args;
} valid_cmds[] = {
  { "declaration",  write_function_declaration, CMD_DECL,   CMDT_SINGLE, 0 },
  { "error",        write_function_error,       CMD_ERROR,  CMDT_ANY,    0 },
  { "nodeclaration",write_function_nodecl,      CMD_NODECL, CMDT_SINGLE, 0 },
  { "typecheck",    write_function_typecheck,   CMD_CHECK,  CMDT_SINGLE, 0 },
  { "return",       write_function_result,      CMD_RETURN, CMDT_EQ,     0 },
//{ "return",       write_function_return,      CMD_RETURN, CMDT_SINGLE, 1 },
  { "return",       write_function_return,      CMD_RETURN, CMDT_0,     1 },
  { "return",       write_function_return,      CMD_RETURN, CMDT_ANY,   1 },
  { "singularcmd",  write_function_singularcmd, CMD_SINGULAR, CMDT_ANY, 1 },
  { NULL,           0, CMD_NONE, CMDT_ANY, 0 }
};

cmd_token checkcmd(
  char *cmdname,
  void (**write_cmd)(moddefv module, procdefv pi, void *arg),
  cmd_type  type,
  int args
  )
{
  int i;
  cmd_token rc = CMD_NONE;
  
  for(i=0; valid_cmds[i].name!=NULL; i++)
    if(strcmp(valid_cmds[i].name, cmdname)==0) {
      rc = CMD_BADSYNTAX;
      if(valid_cmds[i].type == type) {
        *write_cmd = valid_cmds[i].write_cmd;
        return valid_cmds[i].id;
      }
    }
  return rc;
}
  
/*========================================================================*/
struct valid_vars_def {
  char *name;
  var_type type;
  var_token id;
  void (*write_cmd)(moddefv module, var_token type = VAR_NONE,
                    idtyp t, void *arg1 = NULL, void *arg2 = NULL);
} valid_vars[] = {
  { "help",         VAR_STRING,  VAR_HELP,     write_main_variable },
  { "info",         VAR_STRING,  VAR_INFO,     write_main_variable },
  { "package",      VAR_STRING,  VAR_MODULE,   0 },
  { "version",      VAR_STRING,  VAR_VERSION,  write_main_variable },
  { "category",	    VAR_STRING,  VAR_CATEGORY, write_main_variable },
  { NULL,           VAR_UNKNOWN, VAR_NONE, 0 }
};

var_token checkvar(
  char *varname,
  var_type type,
  void (**write_cmd)(moddefv module, var_token type,
                    idtyp t, void *arg1, void *arg2)
  )
{
  int i;
  for(i=0; valid_vars[i].name!=NULL; i++)
    if((strcmp(valid_vars[i].name, varname)==0) &&
       (valid_vars[i].type == type) ) {
      *write_cmd = valid_vars[i].write_cmd;
      return valid_vars[i].id;
    }
  return VAR_NONE;
}

/*========================================================================*/
void write_main_variable(
  moddefv module,
  var_token type,
  idtyp t, 
  void *arg1,
  void *arg2
  )
{
  enter_id(module->fmtfp, t, (char *)arg1, (char *)arg2, yylineno,
           module->filename);
  switch(type) {
    case VAR_INFO:
          module->info = (char *)malloc(strlen((char *)arg2)+1);
	  memset(module->info, '\0', strlen((char *)arg2)+1);
	  memcpy(module->info,(char *)arg2,strlen((char *)arg2));
          break;
    case VAR_VERSION:
          module->version = (char *)malloc(strlen((char *)arg2)+1);
	  memset(module->version, '\0', strlen((char *)arg2)+1);
          memcpy(module->version,(char *)arg2,strlen((char *)arg2));
	  break;
    case VAR_CATEGORY:
          module->category = (char *)malloc(strlen((char *)arg2)+1);
	  memset(module->category, '\0', strlen((char *)arg2)+1);
	  memcpy(module->category,(char *)arg2,strlen((char *)arg2));
	  break;
    default: break;
  }
}
  
/*========================================================================*/
void PrintProc(
  procdefv pi
  )
{
  int i;
  
  printf("%4d proc: %s(", pi->lineno, pi->procname);
  
  for(i=0; i<pi->paramcnt; i++) {
    printf("%s (%s)", (pi->param[i]).name, (pi->param[i]).typname);
    if(i < (pi->paramcnt-1)) printf(",");
  }
  printf(")\n");
  if(pi->return_val.typ!=0)
    printf("\treturn = %s (%s)\n", pi->return_val.name,
           pi->return_val.typname);
  printf("{%s}\n", pi->c_code);
}

/*========================================================================*/
void make_version(char *p, moddefv module)
{
  char ver[10];
  char date[16];
  char libnamebuf[128];
  
  module->major = 0;
  module->minor = 0;
  module->level = 0;
  
  ver[0]='0'; ver[1]='.'; ver[2]='0'; ver[3]='.'; ver[4]='0'; ver[5]='\0';
  date[0]='?'; date[1]='\0';
  sscanf(p,"%*s %*s %10s %16s",ver,date);
  sscanf(ver, "%d.%d.%d", &module->major, &module->minor, &module->level);
  
  sprintf(libnamebuf,"(%s,%s)", ver, date);
  if(strcmp(libnamebuf, "(0.0.0,?)")==0)
  {
    sscanf(p,"%*[^\"]\"%[^\"]\"",libnamebuf);
  }
  module->revision = (char *)malloc(strlen(libnamebuf)+1);
  memset(module->revision, '\0', strlen(libnamebuf)+1);
  memcpy(module->revision, libnamebuf, strlen(libnamebuf));
}

/*========================================================================*/
void make_module_name(char *p, moddefv module)
{
  if(strlen(p)>=1)
  {
    module->targetname = (char *)malloc(strlen(p)+1);
    memset(module->targetname, '\0', strlen(p)+1);
    memcpy(module->targetname,p,strlen(p));
  } else {
    module->targetname = (char *)malloc(strlen(module->name)+1);
    memset(module->targetname, '\0', strlen(module->name)+1);
    memcpy(module->targetname,module->name,strlen(module->name));
  }

}

/*========================================================================*/
void Add2files(
  moddefv module,
  char *name
  )
{
  cfiles cfnew;
  memset((void *)&cfnew, '\0', sizeof(cfiles));

  cfnew.filename = (char *)malloc(strlen(name)+1);
  memset(cfnew.filename, '\0', strlen(name)+1);
  memcpy(cfnew.filename, name, strlen(name));

  if(module->filecnt==0) {
    module->files = (cfilesv)malloc(sizeof(cfiles)+1);
  }
  else {
    module->files = (cfilesv)realloc(module->files,
                                   (module->filecnt+1)*sizeof(cfiles));
  }
  if(module->files == NULL) { printf("ERROR\n"); return; }
  
  memset((void *) &module->files[module->filecnt], '\0', sizeof(cfiles));
  memcpy((void *)(&(module->files[module->filecnt])),
         (void *)&cfnew, sizeof(cfiles));
  (module->filecnt)++;
}

/*========================================================================*/
void init_system_type()
{
  if(strcmp(S_UNAME, "ix86-Linux") == 0) {
    systyp = SYSTYP_LINUX;
  } else if (strcmp(S_UNAME, "HPUX-9")==0) {
    systyp = SYSTYP_HPUX9;
  } else if (strcmp(S_UNAME, "HPUX-10")==0) {
    systyp = SYSTYP_HPUX10;
  }
}

/*========================================================================*/
void PrintProclist(
  moddefv module
  )
{
  logx("PrintProclist()\n");
  int j;
  for(j=0; j<module->proccnt; j++) {
    PrintProc(&(module->procs[j]));
  }
}

/*========================================================================*/
void generate_mod(
  moddefv module,
  int section
  )
{
  procdefv v = NULL;
  cfilesv c_filelist = NULL;
  int proccnt;
  FILE *fp_h;
  char *filename;
  
  switch(section) {
      case 1:
        write_intro(module);
        return;
        
      case 2:
        printf("Writing %s, section %d", filename, section);fflush(stdout);
        break;
  }
  
  sprintf(filename, "%s.h", module->name);
  fp_h = fopen(filename, "w");
  write_header(fp_h, module->name);
  printf("%s  ...", filename);fflush(stdout);
  
  /* building mod_init() */
  
  for(proccnt=0; proccnt<module->proccnt; proccnt++) {
    printf("->%s, %s\n", module->procs[proccnt].procname,
           module->procs[proccnt].funcname);
    fprintf(module->modfp, "  iiAddCproc(\"%s\",\"%s\",%s, mod_%s);\n",
	    module->name, module->procs[proccnt].procname, 
	    module->procs[proccnt].is_static ? "TRUE" : "FALSE",
            module->procs[proccnt].funcname);
  }
  fprintf(module->modfp, "  return 0;\n}\n\n");

  /* building entry-functions */
  for(proccnt=0; proccnt<module->proccnt; proccnt++) {
    //generate_function(&module->procs[proccnt], module->modfp);
    //generate_header(&module->procs[proccnt], fp_h);
  }
  printf("  done.\n");fflush(stdout);
  fclose(module->modfp);
  fclose(fp_h);
}




/*========================================================================*/
void  write_procedure_text(
  moddefv module,
  int lineno
)
{
  int i;
  procdefv pi = &module->procs[module->proccnt-1];

  fprintf(module->fmtfp, "#line %d \"%s\"\n", lineno, module->filename);
  fprintf(module->fmtfp, "%s\n", module->procs[module->proccnt-1].c_code);
  free(module->procs[module->proccnt-1].c_code);
  module->procs[module->proccnt-1].c_code = NULL;

#if 0
  if(pi->paramcnt>0) {
  }
  else {
    switch( pi->return_val.typ) {
        case SELF_CMD:
          fprintf(module->fmtfp, "  return(%s(res));", pi->funcname);
          break;

        case NONE:
          fprintf(module->fmtfp, "  return FALSE;");
          break;
            
        default:
          fprintf(module->fmtfp, "  res->rtyp = %s;\n", 
                  pi->return_val.typname);
          fprintf(module->fmtfp, "  res->data = (void *)%s();\n", 
                  pi->funcname);
          fprintf(module->fmtfp, "  return FALSE;");
    }
  }
#endif
  fprintf(module->fmtfp, "/*\n}\n\n*/");
  
}


/*========================================================================*/
void  mod_write_header(FILE *fp, char *module, char what)
{

  write_header(fp, module);
  fprintf(fp, "#line %d \"%s.cc\"\n", modlineno++, module);
  if(what != 'h') {
    fprintf(fp, "#include <stdlib.h>\n");
    fprintf(fp, "#include <stdio.h>\n");
    fprintf(fp, "#include <string.h>\n");
    fprintf(fp, "#include <ctype.h>\n");
    fprintf(fp, "%s\n", DYNAinclude[systyp]);
    fprintf(fp, "#include <unistd.h>\n");
    fprintf(fp, "#include <sys/stat.h>");
    fprintf(fp, "\n");
    fprintf(fp, "#include <mod2.h>\n");
    fprintf(fp, "#include <tok.h>\n");
    fprintf(fp, "#include <structs.h>\n");
    fprintf(fp, "#include <ipid.h>\n\n");
    fprintf(fp, "#include <locals.h>\n");
    fprintf(fp, "#include <omalloc.h>\n");
    fprintf(fp, "#include \"%s.h\"\n", module);
    modlineno+=8;
    
    fprintf(fp, "#line %d \"%s.cc\"\n", modlineno++, module);
    write_enter_id(fp);
    fprintf(fp, "\n");    modlineno+=1;
    fprintf(fp, "#line %d \"%s.cc\"\n", modlineno++, module);
    write_add_singular_proc(fp);
    write_crccheck(fp);
    fprintf(fp, "\n");
    
    fprintf(fp, "void fill_help_package();\n");
    fprintf(fp, "void fill_example_package();\n");
    modlineno+=3;
  }
  fprintf(fp, "\n");
  
}

/*========================================================================*/
void write_header(FILE *fp, char *module, char *comment)
{
  fprintf(fp, "%s/*\n%s * This was automatically generated by modgen\n",
          comment, comment);
  fprintf(fp, "%s * version %s\n", comment, MOD_GEN_VERSION);
  fprintf(fp, "%s * module %s\n", comment, module);
  fprintf(fp, "%s * Don't edit this file\n%s */\n", comment, comment);
  fprintf(fp, "%s\n", comment);
  fprintf(fp, "%s\n", comment);
  if(strlen(comment)==1)modlineno+=10;
}

/*========================================================================*/
void enter_id(
  FILE *fp, 
  idtyp t,
  char *name,
  char *value,
  int lineno,
  char *file
  )
{
  char tname[32];
  
  if(lineno)
    fprintf(fp, "#line %d \"%s\"\n", lineno, file);
  fprintf(fp, "  enter_id(\"%s\",\"%s\", %s);\n",name, value,
          decl2str(t, tname));
}

/*========================================================================*/
void init_type_conv() 
{
  strcpy(type_conv[NONE], "none");
//  strcpy(type_conv[NONE], "void");
  strcpy(type_conv[INT_CMD], "int");
  strcpy(type_conv[RING_CMD], "ring");
  strcpy(type_conv[QRING_CMD], "ring");
  strcpy(type_conv[POLY_CMD], "poly");
  strcpy(type_conv[NUMBER_CMD], "number");
  strcpy(type_conv[MODUL_CMD], "ideal");
  strcpy(type_conv[VECTOR_CMD], "ideal");
  strcpy(type_conv[IDEAL_CMD], "ideal");
  strcpy(type_conv[MAP_CMD], "map");
  strcpy(type_conv[MATRIX_CMD], "matrix");
  strcpy(type_conv[STRING_CMD], "char *");
  strcpy(type_conv[INTMAT_CMD], "intvec *");
  strcpy(type_conv[INTVEC_CMD], "intvec *");
  strcpy(type_conv[LIST_CMD], "lists");
  strcpy(type_conv[LINK_CMD], "si_link");
  strcpy(type_conv[PACKAGE_CMD], "package");
  strcpy(type_conv[PROC_CMD], "procinfo *");
/*
  strcpy(type_conv[], "");
  strcpy(type_conv[], "");
  strcpy(type_conv[], "");
  strcpy(type_conv[], "");
  strcpy(type_conv[], "");
  strcpy(type_conv[], "");
  strcpy(type_conv[], "");
  printf("[%d] %s\n", INT_CMD, type_conv[INT_CMD]);
  printf("[%d] %s\n", MODUL_CMD, type_conv[MODUL_CMD]);
  printf("[%d] %s\n", STRING_CMD, type_conv[STRING_CMD]);
  printf("[%d] %s\n", LINK_CMD, type_conv[LINK_CMD]);
  printf("[%d] %s\n", PACKAGE_CMD, type_conv[PACKAGE_CMD]);
  / **/
}

/*========================================================================*/
