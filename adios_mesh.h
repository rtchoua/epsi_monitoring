struct PairStrct {
    char * name;
    char * value;
    struct PairStrct * next;
};
typedef struct PairStrct PairStrct;

struct schema_mesh_item_struct
{
    double rank;
    struct adios_var_struct * var;
};

struct schema_mesh_item_list_struct
{
    struct schema_mesh_item_struct item;
    struct schema_mesh_item_list_struct * next;
};

struct schema_mesh_var_list_struct
{
    struct adios_var_struct * var;
    struct schema_mesh_var_list_struct * next;
};

struct schema_mesh_cell_list_struct
{
    enum ADIOS_FLAG cells_uniform;
    struct schema_mesh_item_struct count;
    struct adios_var_struct * data;
    struct schema_mesh_item_struct type;
};

struct schema_mesh_cell_list_list_struct
{
    struct schema_mesh_cell_list_struct cell_list;
    struct schema_mesh_cell_list_list_struct * next;
};

enum SCHEMA_FLAG {
         SCHEMA_FLAG_UNKNOWN = 0
        ,SCHEMA_FLAG_YES = 1
        ,SCHEMA_FLAG_NO = 2
};

enum SCHEMA_MESH_TIME_VARYING
{
     SCHEMA_MESH_TIME_VARY_YES  = 1
    ,SCHEMA_MESH_TIME_VARY_NO   = 0
};

enum SCHEMA_MESH_FILE
{
     SCHEMA_MESH_FILE_YES  = 1
    ,SCHEMA_MESH_FILE_NO   = 0
}
;
enum SCHEMA_MESH_GROUP
{
     SCHEMA_MESH_GROUP_YES  = 1
    ,SCHEMA_MESH_GROUP_NO   = 0
};

enum SCHEMA_MESH_TYPE
{
     SCHEMA_MESH_UNIFORM      = 1
    ,SCHEMA_MESH_STRUCTURED   = 2
    ,SCHEMA_MESH_RECTILINEAR  = 3
    ,SCHEMA_MESH_UNSTRUCTURED = 4
};

// ADIOS Schema: supported cell types
enum SCHEMA_CELL_TYPE
{
     SCHEMA_CELL_PT         = 1
    ,SCHEMA_CELL_LINE       = 2
    ,SCHEMA_CELL_TRI        = 3
    ,SCHEMA_CELL_QUAD       = 4
    ,SCHEMA_CELL_HEX        = 5
    ,SCHEMA_CELL_PRI        = 6
    ,SCHEMA_CELL_TET        = 7
    ,SCHEMA_CELL_PYR        = 8
};

struct schema_meshstruct;

////////////////////////
// mesh support data structures
////////////////////////
struct schema_meshitem_struct
{
    double rank;
    struct adios_var_struct * var;
};

struct schema_meshitem_list_struct
{
    struct schema_meshitem_struct item;
    struct schema_meshitem_list_struct * next;
};

struct schema_meshvar_list_struct
{
    struct adios_var_struct * var;
    struct schema_meshvar_list_struct * next;
};

struct schema_meshcell_list_struct
{
    enum SCHEMA_FLAG cells_uniform;
    struct schema_meshitem_struct count;
    struct adios_var_struct * data;
    struct schema_meshitem_struct type;
};
struct schema_meshcell_list_list_struct
{
    struct schema_meshcell_list_struct cell_list;
    struct schema_meshcell_list_list_struct * next;
};

//////////////////////////////////////////////////////////////
// Main mesh structs
//////////////////////////////////////////////////////////////
struct schema_mesh_struct
{
    // ADIOS Schema: adding mesh names 
    // Groups can have multiple meshes
    char * name;
    char * filename;
    char * groupname;
    int time_series_format;
    enum SCHEMA_FLAG time_varying;
    enum SCHEMA_MESH_TYPE type;
    enum SCHEMA_MESH_FILE file;
    enum SCHEMA_MESH_GROUP group;
    union
    {
        struct schema_mesh_uniform_struct * uniform;
        struct schema_mesh_rectilinear_struct * rectilinear;
        struct schema_mesh_structured_struct * structured;
        struct schema_mesh_unstructured_struct * unstructured;
    };
    struct schema_meshstruct * next;
};

struct schema_mesh_uniform_struct
{
    struct schema_meshitem_list_struct * dimensions;
    struct schema_meshitem_list_struct * origin;
    struct schema_meshitem_list_struct * spacing;
    // ADIOS Schema: adding option to provide origin and maximum
    // instead restricting users to origin and spacing
    struct schema_meshitem_list_struct * maximum;
};

struct schema_mesh_rectilinear_struct
{
    enum SCHEMA_FLAG coordinates_single_var;
    struct schema_meshitem_list_struct * dimensions;
    struct schema_meshvar_list_struct * coordinates;
};

struct schema_mesh_structured_struct
{
    enum SCHEMA_FLAG points_single_var;
    struct schema_meshitem_struct * nspace;
    struct schema_meshitem_list_struct * dimensions;
    struct schema_meshvar_list_struct * points;
};

struct schema_mesh_unstructured_struct
{
    // ADIOS Schema: adding single/multi points option
    // adding nspace to allow 2D mesh in 3D for example,
    // finally adding the concept of cellset/cellsetcount
    enum SCHEMA_FLAG points_single_var;
    struct schema_meshitem_struct * nspace;
    struct schema_meshvar_list_struct * points;
    struct schema_meshitem_struct * points_count;
    struct schema_meshcell_list_list_struct * cell_list;
    struct schema_meshitem_struct * cell_set_count;
};


