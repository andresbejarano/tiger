#include "message.h"

Message::Message() 
{
    set(MSG_DEFAULT);
}

Message::Message(MESSAGE_TYPE type)
{
    set(type);
}

Message::~Message() 
{
}

const char * Message::getText() const
{
    switch (m_type)
    {
        case MSG_UNKNOWN:               return STR_UNKNOWN;             break;
        case MSG_BLOCKS_MODIFIED:       return STR_BLOCKS_MODIFIED;     break;
        case MSG_NO_BLOCKS:             return STR_NO_BLOCKS;           break;
        case MSG_NO_EDGE_DIRECTIONS:	return STR_NO_EDGE_DIRECTIONS;  break;
        case MSG_NO_EVEN_SIDED_TILES:   return STR_NO_EVEN_SIDED_TILES; break;
        case MSG_NO_INTERFACES:         return STR_NO_INTERFACES;       break;
        case MSG_NO_TESSELLATION:       return STR_NO_TESSELLATION;     break;
        case MSG_NO_TILE_CENTERS:       return STR_NO_TILE_CENTERS;     break;
        case MSG_OK:                    return STR_OK;                  break;
        case MSG_DEFAULT:               return STR_DEFAULT;             break;
        default:                        return STR_UNKNOWN;
    }
}

Message::MESSAGE_TYPE Message::getType() const 
{
    return m_type;
}

void Message::set(MESSAGE_TYPE type)
{
    m_type = (type < MSG_UNKNOWN || type > MSG_DEFAULT) ? MSG_UNKNOWN : type;
}
