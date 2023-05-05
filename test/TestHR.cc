#include <StableMarriage.h>
#include <Popular.h>
#include <Utils.h>
#include "TestDefs.h"

TEST_CASE("StableMarriage hrlq_m6", "[matching_SM_HR]") {
    auto G = read_graph(get_filepath(get_resources_dir(), "/hrlq_m6.txt"));
    auto r1 = get_vertex_by_id(G, "r1");
    auto r2 = get_vertex_by_id(G, "r2");
    auto r3 = get_vertex_by_id(G, "r3");
    auto r4 = get_vertex_by_id(G, "r4");

    auto h1 = get_vertex_by_id(G, "h1");
    auto h2 = get_vertex_by_id(G, "h2");
    auto h3 = get_vertex_by_id(G, "h3");

    SECTION("residents proposing") {
        StableMarriage sm(G);
        auto M = sm.compute_matching();

        REQUIRE(M.size() == 3);

        SECTION("size of partner list") {
            REQUIRE(M.number_of_partners(r1) == 1);
            REQUIRE(M.number_of_partners(r2) == 1);
            REQUIRE(M.number_of_partners(r3) == 1);
            REQUIRE(M.number_of_partners(h1) == 2);
            REQUIRE(M.number_of_partners(h2) == 1);
        }

        SECTION("actual partners") {
            REQUIRE(M.get_partner(r1) == h1);
            REQUIRE(M.get_partner(r2) == h1);
            REQUIRE(M.get_partner(r3) == h2);
            REQUIRE(M.has_partners(h1, {r1, r2}));
            REQUIRE(M.has_partners(h2, {r3}));
        }
    }

    SECTION("hospitals proposing") {
        StableMarriage sm(G, false);
        auto M = sm.compute_matching();

        REQUIRE(M.size() == 3);

        SECTION("size of partner list") {
            REQUIRE(M.number_of_partners(r1) == 1);
            REQUIRE(M.number_of_partners(r2) == 1);
            REQUIRE(M.number_of_partners(r3) == 1);
            REQUIRE(M.number_of_partners(h1) == 2);
            REQUIRE(M.number_of_partners(h2) == 1);
        }

        SECTION("actual partners") {
            REQUIRE(M.get_partner(r1) == h1);
            REQUIRE(M.get_partner(r2) == h1);
            REQUIRE(M.get_partner(r3) == h2);
            REQUIRE(M.has_partners(h1, {r1, r2}));
            REQUIRE(M.has_partners(h2, {r3}));
        }
    }
}

TEST_CASE("MaxCardPopular hrlq_m6", "[matching_MP_HR]") {
    auto G = read_graph(get_filepath(get_resources_dir(), "/hrlq_m6.txt"));
    auto r1 = get_vertex_by_id(G, "r1");
    auto r2 = get_vertex_by_id(G, "r2");
    auto r3 = get_vertex_by_id(G, "r3");
    auto r4 = get_vertex_by_id(G, "r4");

    auto h1 = get_vertex_by_id(G, "h1");
    auto h2 = get_vertex_by_id(G, "h2");
    auto h3 = get_vertex_by_id(G, "h3");

    SECTION("residents proposing") {
        MaxCardPopular mp(G);
        auto M = mp.compute_matching();

        REQUIRE(M.verify(G));
        REQUIRE(M.size() == 4);

        SECTION("size of partner list") {
            REQUIRE(M.number_of_partners(r1) == 1);
            REQUIRE(M.number_of_partners(r2) == 1);
            REQUIRE(M.number_of_partners(r3) == 1);
            REQUIRE(M.number_of_partners(h1) == 2);
            REQUIRE(M.number_of_partners(h2) == 1);
            REQUIRE(M.number_of_partners(h3) == 1);
        }

        SECTION("actual partners") {
            REQUIRE(M.get_partner(r1) == h2);
            REQUIRE(M.get_partner(r2) == h3);
            REQUIRE(M.get_partner(r3) == h1);
            REQUIRE(M.get_partner(r4) == h1);
            REQUIRE(M.has_partners(h1, {r3, r4}));
            REQUIRE(M.has_partners(h2, {r1}));
            REQUIRE(M.has_partners(h3, {r2}));
        }
    }

    SECTION("hospitals proposing") {
        MaxCardPopular mp(G, false);
        auto M = mp.compute_matching();

        REQUIRE(M.verify(G));
        REQUIRE(M.size() == 4);

        SECTION("size of partner list") {
            REQUIRE(M.number_of_partners(r1) == 1);
            REQUIRE(M.number_of_partners(r2) == 1);
            REQUIRE(M.number_of_partners(r3) == 1);
            REQUIRE(M.number_of_partners(r4) == 1);
            REQUIRE(M.number_of_partners(h1) == 2);
            REQUIRE(M.number_of_partners(h2) == 1);
            REQUIRE(M.number_of_partners(h3) == 1);
        }

        SECTION("actual partners") {
            REQUIRE(M.get_partner(r1) == h3);
            REQUIRE(M.get_partner(r2) == h1);
            REQUIRE(M.get_partner(r3) == h1);
            REQUIRE(M.get_partner(r4) == h2);
            REQUIRE(M.has_partners(h1, {r2, r3}));
            REQUIRE(M.has_partners(h2, {r4}));
            REQUIRE(M.has_partners(h3, {r1}));
        }
    }
}

TEST_CASE("MaxCardPopular P_mat_level_up", "[matching_MP_HR_level_up]") {
    auto G = read_graph(get_filepath(get_resources_dir(), "/P_mat_level_up.txt"));
    auto r1 = get_vertex_by_id(G, "r1");
    auto r4 = get_vertex_by_id(G, "r4");
    auto r5 = get_vertex_by_id(G, "r5");
    auto r7 = get_vertex_by_id(G, "r7");
    auto r8 = get_vertex_by_id(G, "r8");
    auto r10 = get_vertex_by_id(G, "r10");

    auto h1 = get_vertex_by_id(G, "h1");
    auto h2 = get_vertex_by_id(G, "h2");
    auto h3 = get_vertex_by_id(G, "h3");

    SECTION("residents proposing") {
        MaxCardPopular mp(G);
        auto M = mp.compute_matching();

        REQUIRE(M.verify(G));
        REQUIRE(M.size() == 6);

        SECTION("size of partner list") {
            REQUIRE(M.number_of_partners(r1) == 1);
            REQUIRE(M.number_of_partners(r4) == 1);
            REQUIRE(M.number_of_partners(r5) == 1);
            REQUIRE(M.number_of_partners(r7) == 1);
            REQUIRE(M.number_of_partners(r8) == 1);
            REQUIRE(M.number_of_partners(r10) == 1);
            REQUIRE(M.number_of_partners(h1) == 2);
            REQUIRE(M.number_of_partners(h2) == 2);
            REQUIRE(M.number_of_partners(h3) == 2);
        }

        SECTION("actual partners") {
            REQUIRE(M.get_partner(r1) == h3);
            REQUIRE(M.get_partner(r4) == h2);
            REQUIRE(M.get_partner(r5) == h1);
            REQUIRE(M.get_partner(r7) == h3);
            REQUIRE(M.get_partner(r8) == h1);
            REQUIRE(M.get_partner(r10) == h2);
            REQUIRE(M.has_partners(h1, {r5, r8}));
            REQUIRE(M.has_partners(h2, {r4, r10}));
            REQUIRE(M.has_partners(h3, {r1, r7}));
        }
    }

    SECTION("hospitals proposing") {
        MaxCardPopular mp(G, false);
        auto M = mp.compute_matching();

        REQUIRE(M.verify(G));
        REQUIRE(M.size() == 6);

        SECTION("size of partner list") {
            REQUIRE(M.number_of_partners(r1) == 1);
            REQUIRE(M.number_of_partners(r4) == 1);
            REQUIRE(M.number_of_partners(r5) == 1);
            REQUIRE(M.number_of_partners(r7) == 1);
            REQUIRE(M.number_of_partners(r8) == 1);
            REQUIRE(M.number_of_partners(r10) == 1);
            REQUIRE(M.number_of_partners(h1) == 2);
            REQUIRE(M.number_of_partners(h2) == 2);
            REQUIRE(M.number_of_partners(h3) == 2);
        }

        SECTION("actual partners") {
            REQUIRE(M.get_partner(r1) == h3);
            REQUIRE(M.get_partner(r4) == h2);
            REQUIRE(M.get_partner(r5) == h1);
            REQUIRE(M.get_partner(r7) == h3);
            REQUIRE(M.get_partner(r8) == h1);
            REQUIRE(M.get_partner(r10) == h2);
            REQUIRE(M.has_partners(h1, {r5, r8}));
            REQUIRE(M.has_partners(h2, {r4, r10}));
            REQUIRE(M.has_partners(h3, {r1, r7}));
        }
    }
}
