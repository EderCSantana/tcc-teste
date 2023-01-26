#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

std::vector<int> DDP(
    int D,
    const std::vector<int>& d )
{
    std::vector<int> P;          // calculated discretization points

    std::vector<int> c( D+1 );

    for( int i = 0; i < (int)d.size(); i++ )
    {
        for( int j = d[i]; j <= D; j++ )
        {
            if( c[ j ] < c[ j - d[i] ] + d[ i ] )
                c[ j ] = c[ j - d[i] ] + d[ i ];
        }
    }
    for( int j = 1; j <= D; j++ )
    {
        if( c[ j ] == j )
            P.push_back( j );
    }
    return P;
}
enum class eCut
{
    nil,
    vert,
    depth,
    horz,
};
std::string text( eCut cut );

typedef std::vector<std::vector<std::vector<std::vector<int> > > > v4d_t;

/// A 3D object with dimension and location
class cSpace
{
public:
    int myLength, myWidth, myHeight;    // dimensions
    int myLocL, myLocW, myLocH;


    cSpace( int L, int W, int H )
        : myLength( L ), myWidth( W ), myHeight{ H }
    {

    }
    bool fit( const cSpace& squeeze ) const
    {
        return ( squeeze.myLength <= myLength &&
                 squeeze.myLength <= myLength &&
                 squeeze.myHeight <= myHeight );
    }

    int size_horiz()
    {
        return myLength * myWidth;
    }
};

// An item that packs into a bin, specialization of cSpace
class cItem : public cSpace
{
public:
    int demand;
    bool myPacked;

    cItem( int L, int W, int H, int D )
        : cSpace( L, W, H )
        , demand( D )
        , myPacked( false )
    {

    }
    void pack( const cItem& space )
    {
        myPacked = true;
        myLocL = space.myLocL;
        myLocW = space.myLocW;
        myLocH = space.myLocH;
    }
};
typedef std::vector< cItem > itemv_t;
std::vector<int> DDP(
    int D,
    const std::vector<int>& d );
std::vector<int> RRP(
    int D,
    const std::vector<int>& d );
bool RK2FFG(
    itemv_t& level,
    const std::vector<int>& bin );
struct sInstance
{
    std::vector<int> bin;
    std::vector<int> l;
    std::vector<int> w;
    std::vector<int> h;
    std::vector<int> demand;        // number of items required to pack
    std::vector<int> item_values;
    std::string myName;
    int stageCount;                 // the number of stages
    void read( const std::string& fname );
    std::string text();
};

struct sPattern
{
    sInstance instance;
    std::vector<int> l_raster, w_raster, h_raster;          /// points where cuts might be made
    v4d_t value;
    v4d_t position;
    v4d_t item;
    v4d_t direction;                         /// cut direction

    std::string text() const;
    std::string textCuts( int stage, eCut cut ) const;
};

sPattern DP3UK (
    sInstance& problem );
sPattern DPS3UK (
    sInstance& problem );
sPattern H3CS (
    sInstance& problem );
std::vector< int >
ParseSpaceDelimited(
    const std::string& l )
{
    std::vector< int > output;
    std::stringstream sst(l);
    std::string a;
    while( getline( sst, a, ' ' ) )
        output.push_back( atoi(a.c_str()));
    return output;
}

std::string text( eCut cut )
{
    switch( cut )
    {
    case eCut::vert:
        return "Vertical";
    case eCut::depth:
        return "Depth";
    case eCut::horz:
        return "Horizontal";
    default:
        return "nil";
    }
}

void sInstance::read( const std::string& fname )
{
    std::ifstream f( fname );
    if( ! f.is_open() )
        throw std::runtime_error("Cannot read instance file " + fname );

    l.clear();
    w.clear();
    h.clear();

    int box_type_count;
    int lcount = 0;
    std::string line;
    while( std::getline( f, line ) )
    {
        if( lcount == 0 )
        {
            box_type_count = atoi( line.c_str() );
            lcount++;
            continue;
        }
        if( lcount == 1 )
        {
            bin = ParseSpaceDelimited( line );
            lcount++;
            continue;          
        }
        if( lcount < 2 + box_type_count )
        {
            // parse an item
            auto lv = ParseSpaceDelimited( line );
            l.push_back( lv[0] );
            w.push_back( lv[1] );
            h.push_back( lv[2] );
            item_values.push_back( lv[0]*lv[1]*lv[2] );
            if( lv.size() == 4 )
                demand.push_back(lv[3]);

            lcount++;
            continue;
        }
        myName = line;
    }
    std::cout << "read problem " << myName << "\n";
    
}

std::string sInstance::text()
{
    std::stringstream ss;
    for( int k = 0; k < (int)l.size(); k++ )
        ss << l[k] <<" "<< w[k] <<" " << h[k] << "\n";
    return ss.str();
}


std::string sPattern::text() const
{
    std::stringstream ss;

    int lCount = l_raster.size();
    int wCount = w_raster.size();
    int hCount = h_raster.size();

    int totalCuts = 0;
    if( instance.stageCount )
    {
        // staged solution
        for( int kstage = 1; kstage < instance.stageCount; kstage++ )
        {
            int stageCuts = 0;
            for( int il = 0; il < lCount; il++ )
            {
                for( int iw = 0; iw < wCount; iw++ )
                {
                    for( int ih = 0; ih < hCount; ih++ )
                    {
                        if( direction[kstage][il][iw][ih]  )
                            stageCuts++;
                    }
                }
            }
            if( instance.stageCount )
            {
                ss << textCuts( kstage, eCut::vert);
                ss << textCuts( kstage, eCut::depth);
                ss << textCuts( kstage, eCut::horz);
                ss << instance.bin[kstage] << "\n";
            }
            totalCuts += stageCuts;
        }
    }
    else
    {
        // unstaged
        for( int il = 0; il < lCount; il++ )
        {
            for( int iw = 0; iw < wCount; iw++ )
            {
                for( int ih = 0; ih < hCount; ih++ )
                {
                    if( direction[0][il][iw][ih]  )
                        totalCuts++;
                }
            }
        }
        ss << textCuts( 0, eCut::vert);
        ss << textCuts( 0, eCut::depth);
        ss << textCuts( 0, eCut::horz);
        ss << "\n";
    }

    if( ! totalCuts )
    {
        ss << "\nno cuts found\n"
           << "Probably means that items were too big for bin\n";
        return ss.str();
    }

    int itemCount = 0;
    int totalValue = 0;
    for( int il = 0; il < lCount; il++ )
    {
        for( int iw = 0; iw < wCount; iw++ )
        {
            for( int ih = 0; ih < hCount; ih++ )
            {
                itemCount++;
                totalValue += instance.item_values[ item[0][il][iw][ih] ];
            }
        }
    }
    
    return ss.str();
}
std::string sPattern::textCuts( int stage, eCut cut ) const
{
    std::stringstream ss;
    std::set<int> cutset;
    for( int il = 0; il < (int)l_raster.size(); il++ )
    {
        for( int iw = 0; iw < (int)w_raster.size(); iw++ )
        {
            for( int ih = 0; ih < (int)h_raster.size(); ih++ )
            {
                if( direction[stage][il][iw][ih] == (int)cut )
                {
                    cutset.insert( position[stage][il][iw][ih] );
                }
            }
        }
    }
    if( cutset.size() )
    {
        if( stage )
            ss << "Stage " << stage << " ";

        ss << cutset.size() << " " << ::text( cut ) << " cuts at ";
        for( int c : cutset )
        {
            ss << c << " ";
        }
        ss << "\n";
    }
    return ss.str();
}
sPattern DPS3UK (
    sInstance& problem )
{
    int L = problem.bin[0];
    int W = problem.bin[1];
    int H = problem.bin[2];
    std::vector<int> l = problem.l;
    std::vector<int> w = problem.w;
    std::vector<int> h = problem.h;
    std::vector<int> v = problem.item_values;
    int k = problem.stageCount;

    eCut previous;

    auto Phat = DDP( L, l );
    auto Qhat = DDP( W, w );
    auto Rhat = DDP( H, h );
    int m = Phat.size();
    int s = Qhat.size();
    int u = Rhat.size();
    std::vector<int>                                G1(u,0);
    std::vector<std::vector<int>  >                 G2(s,G1);
    std::vector<std::vector<std::vector<int> > >    G3(m,G2);
    v4d_t    G(k,G3);        // the best value for boxes in this bin
    v4d_t    item(m,G3);     // the item at this location
    v4d_t    guil(m,G3);     // the cut orientation
    v4d_t    pos(m,G3);      // the cut position
    for( int ilength = 0; ilength < m; ilength++ )
    {
        for( int iwidth = 0; iwidth < s; iwidth++ )
        {
            for( int iheight = 0; iheight < u; iheight++ )
            {
                // loop over the items
                for( int d = 0; d < (int)v.size(); d++ )
                {
                    if( l[d] <= Phat[ilength] &&
                            w[d] <= Qhat[iwidth] &&
                            h[d] <= Rhat[iheight] )
                    {
                        // the item is small enough to fit

                        // does it have greater value than previous fitted items?
                        if( v[d] > G[0][ilength][iwidth][iheight] )
                        {
                            G[0][ilength][iwidth][iheight] = v[d];
                            item[0][ilength][iwidth][iheight] = d;
                            guil[0][ilength][iwidth][iheight] = (int) eCut::nil;
                        }
                    }
                }

            }
        }
    }
    //pseudocode line 2.9
    switch( k % 3 )
    {
    case 1:
        previous = eCut::vert;
        break;
    case 2:
        previous = eCut::depth;
        break;
    default:
        previous = eCut::horz;
        break;
    }

    // loop over stages
    for( int b = 1; b < k; b++ )
    {
        for( int ilength = 0; ilength < m; ilength++ )
        {
            for( int iwidth = 0; iwidth < s; iwidth++ )
            {
                for( int iheight = 0; iheight < u; iheight++ )
                {
                    G[b][ilength][iwidth][iheight] = G[b-1][ilength][iwidth][iheight];
                    guil[b][ilength][iwidth][iheight] = (int)eCut::nil;
                    if( previous == eCut::depth )
                    {
                        // avoid generating symmetric patterns by considering, in each direction,
                        // r-points up to half of the size of the respective bin
                        int nn = -1;
                        for( int d = 0; d <= ilength; d++ )
                        {
                            if( Phat[d] <= Phat[ilength] / 2 )
                            {
                                nn = d;
                            }
                        }
                        for( int x = 0; x <= nn; x++ )
                        {
                            int t = 0;
                            for( int d = 0; d < m; d++ )
                            {
                                if( Phat[d] <= Phat[ilength] - Phat[x] )
                                    t = d;
                            }
                            // Does vertical cut at Phat[x] improve value of solution
                            if( G[b][ilength][iwidth][iheight] < G[b][x][iwidth][iheight]+G[b][t][iwidth][iheight] )
                            {


                                G[b][ilength][iwidth][iheight] = G[b][x][iwidth][iheight]+G[b][t][iwidth][iheight];
                                pos[b][ilength][iwidth][iheight] = Phat[t];
                                guil[b][ilength][iwidth][iheight] = (int) eCut::vert;  // Vertical cut; parallel to yz-plane

                            }
                        }
                    }

                    else if ( previous == eCut::vert )
                    {
                        int nn = -1;
                        for( int d = 0; d <= iheight; d++ )
                        {
                            if( Rhat[d] <= Rhat[iheight] / 2 )
                                nn = d;
                        }
                        for( int z = 0; z <= nn; z++ )
                        {
                            int t = 0;
                            for( int d = 0; d < u; d++ )
                            {
                                if( Rhat[d] <= Rhat[iheight] - Rhat[z] )
                                    t = d;
                            }
                            if( G[b][ilength][iwidth][iheight] < G[b][ilength][iwidth][z]+G[b][ilength][iwidth][t] )
                            {


                                G[b][ilength][iwidth][iheight] = G[b][ilength][iwidth][z]+G[b][ilength][iwidth][t];
                                pos[b][ilength][iwidth][iheight] = Rhat[t];
                                guil[b][ilength][iwidth][iheight] = (int) eCut::horz;  // Horizontal cut; parallel to xy plane
                            }
                        }
                    }
                    else
                    {
                        int nn = -1;
                        for( int d = 0; d <= iwidth; d++ )
                        {
                            if( Qhat[d] <= Qhat[iwidth] / 2 )
                                nn = d;
                        }
                        for( int y = 0; y <= nn; y++ )
                        {
                            int t = 0;
                            for( int d = 0; d < s; d++ )
                            {
                                if( Qhat[d] <= Qhat[iwidth] - Qhat[y] )
                                    t = d;
                            }
                            if( G[b][ilength][iwidth][iheight] < G[b][ilength][y][iheight]+G[b][ilength][t][iheight] )
                            {
                                G[b][ilength][iwidth][iheight] = G[b][ilength][y][iheight]+G[b][ilength][y][iheight];
                                pos[b][ilength][iwidth][iheight] = Qhat[t];
                                guil[b][ilength][iwidth][iheight] = (int) eCut::depth;  // Depth cut �vertical; parallel to xy plane
                            }
                        }
                    }
                }
            }
        }
        // next stage uses next cut orientation
        switch( previous )
        {
        case eCut::vert:
            previous = eCut::depth;
            break;
        case eCut::depth:
            previous = eCut::horz;
            break;
        case eCut::horz:
            previous = eCut::vert;
            break;
        default:
            throw std::runtime_error("DPS3UK Bad stage value");
        }

    }   // end loop over stages
    sPattern P;
    P.instance = problem;

    P.value = G ;
    P.position = pos;
    P.direction = guil;
    P.item =item;
    P.l_raster = Phat;
    P.w_raster = Qhat;
    P.h_raster = Rhat;

    return P;
}

int main( int argc, char* argv[] )
{
    std::cout << "DP3SUK\n";

    if( argc != 2 )
    {
        std::cout << "Usage: DP3SUKInstance data/thpack/thpack1\n";
        exit(1);
    }

    sInstance problem;
    problem.read( argv[1] );
    problem.stageCount = 3;
    auto P = DPS3UK( problem );


    std::cout << P.text() << "\n";

    return 0;
}
