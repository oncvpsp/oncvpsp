set(src_wxml
    xmlf90_wxml.f90
    m_wxml_buffer.f90
    m_wxml_core.f90
    m_wxml_dictionary.f90
    m_wxml_elstack.f90
    m_wxml_text.f90
)

###########################################################
# xmlf90-wxml library
###########################################################
add_library(wxml STATIC ${src_wxml})
