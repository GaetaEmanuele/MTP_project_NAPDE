function [tents,nodes,tentDependencies] = create_tents(Dati)
    x0 = Dati.domain(1);
    xN = Dati.domain(2);
    h = Dati.h;
    N = (xN-x0)/h;
    T = Dati.T;
    x = linspace(x0,xN,N);
    c = Dati.c;
    [nodes,n_even,n_odd] = C_nodes(x,T,c,N);
    x_dim = size(nodes,2);
    max_it = size(nodes,1)-2;
    t = 2;
    x_cord = [nodes(1,1),nodes(1,2)];
    height = [nodes(t,1),nodes(t,2),nodes(t+1,1)];
    Tent = StructTent(x_cord, height);
    tents = [];
    tents(:,:,1) = Tent;
    while(t <= max_it-1)
        x_cord = [nodes(1,1),nodes(1,2)];
        if rem(t,2)==0 && t>2
            for i=1:2:(x_dim-3)
                x_cord = [nodes(1,i),nodes(1,i+1),nodes(1,i+2)];
                height = [nodes(t,i),nodes(t-1,i+1),nodes(t,i+2),nodes(t+1,i+1)];
                Tent = StructTent(x_cord, height);
                tents(:,:,end+1) = Tent;
            end
        elseif rem(t,2)>0 && t>2
            for i=2:2:(x_dim-2)
                if(nodes(t-1,i)==0)
                    x_cord = [nodes(1,i),nodes(1,i+1),nodes(1,i+2)];
                    height = [nodes(t,i),nodes(t+1,i+1),nodes(t,i+2)];
                    Tent = StructTent(x_cord, height);
                    tents(:,:,end+1) = Tent;
                else
                    x_cord = [nodes(1,i),nodes(1,i+1),nodes(1,i+2)];
                    height = [nodes(t,i),nodes(t-1,i+1),nodes(t,i+2),nodes(t+1,i+1)];
                    Tent = StructTent(x_cord, height);
                    tents(:,:,end+1) = Tent;
                end
            end
        end
        if t>2 && rem(t,2)==0
           x_cord = [nodes(1,1),nodes(1,2)];
           height = [nodes(t-1,1),nodes(t,2),nodes(t+1,1)];
           Tent = StructTent(x_cord, height);
           tents(:,:,end+1) = Tent;
        end
        if t>2 && rem(t,2)>0
           x_cord = [nodes(1,end-1),nodes(1,end)];
           height = [nodes(t,end-1),nodes(t-1,end),nodes(t+1,end)];
           Tent = StructTent(x_cord, height);
           tents(:,:,end+1) = Tent;
        end
        t = t+1;
    end
    x_cord = [nodes(1,1),nodes(1,2)];
    height = [nodes(t,1),nodes(t-1,1),nodes(t,2)];
    Tent = StructTent(x_cord, height);
    tents(:,:,end+1) = Tent;
    for i=2:2:(x_dim-2)
        x_cord = [nodes(1,i),nodes(1,i+1),nodes(1,i+2)];
        height = [nodes(t,i),nodes(t-1,i+1),nodes(t,i+2)];
        Tent = StructTent(x_cord, height);
        tents(:,:,end+1) = Tent;
    end
    tentDependencies = findTentDependencies(tents,n_even,n_odd);
end