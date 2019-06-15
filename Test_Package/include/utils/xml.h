#ifndef __has_utils_xml_h
#define __has_utils_xml_h

#ifdef ENABLE_XML2

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <string>

namespace utils
{

xmlDoc*		xml_open(const char *filename);
xmlNode*	xml_find_element_doc (xmlDoc *doc, const char *element);
xmlNode*	xml_find_element_node (xmlNode *node, const char *element);
int		xml_close (xmlDoc *doc);
void		xml_print_element_names(xmlNode * a_node, int depth = 0);
int		xml_parse_property (const char *value, xmlNode *node, std::string& str);
int		xml_parse_element_content (xmlNode *node, std::string &str);

} // utils

#endif

#endif
