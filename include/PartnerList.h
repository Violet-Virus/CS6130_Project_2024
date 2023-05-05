#ifndef PARTNER_LIST_H
#define PARTNER_LIST_H

#include "TDefs.h"
#include "Partner.h"
#include <list>
#include <memory>

class PartnerList {
public:
    typedef typename std::list<Partner> ContainerType;
    typedef ContainerType::iterator Iterator;
    typedef ContainerType::const_iterator ConstIterator;
    typedef ContainerType::size_type SizeType;

private:
    ContainerType partners_;

public:
    PartnerList() = default;
    virtual ~PartnerList();

    bool empty() const;
    SizeType size() const;

    ConstIterator cbegin() const noexcept;
    ConstIterator cend() const noexcept;

    Iterator begin() noexcept;
    Iterator end() noexcept;

    ConstIterator begin() const noexcept;
    ConstIterator end() const noexcept;

    /// is this vertex in the list
    ConstIterator find(VertexPtr v) const;

    /// add a vertex to the list of matched partners
    void add_partner(VertexPtr partner, RankType rank, int level);

    // Get the partner from this list. Must have 0 or 1 partners.
    VertexPtr get_partner() const;

    ConstIterator find_least_preferred() const;
    /// return details for the worst partner matched to this vertex
    Partner get_least_preferred() const;

    /// remove this partner from the list
    void remove(VertexPtr v);

    /// remove the least preferred among the current partners
    void remove_least_preferred();

/*
    friend std::ostream& operator<<(std::ostream& out, PartnerList& pl);
    friend std::ostream& operator<<(std::ostream& out, PartnerList* pl);
*/
};

#endif
