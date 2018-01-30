%% Graham's scan
% 14. \textbf{return} L
function L = convex_hull(P)
% 1. Sort the points by x?coordinate, resulting in a sequence
% p1,..., pn.
Ps = sortrows(P,2); % sort points by y?value
Ps = sortrows(Ps,1); % sort points by x?value
% X. Remove duplicate points
for c = length(Ps):-1:2
if isequal(Ps(c,:),Ps(c-1,:)), Ps(c,:) = []; 
end
end
% 2. Put the points p1 and p2 in a list Lupper, with
% p1 as the first point.
Lu = Ps(1:2,:);
% 3. \textbf{for} ii ¡û 3 \textbf{to} n
for ii = 3:length(Ps)
% 4. \textbf{do} Append pi to Lupper.
Lu = [Lu;Ps(ii,:)];
% 5.1 \textbf{while} Lupper contains more than 2 points
while length(Lu) > 2
p1 = [Lu(end-2,:),0];
p2 = [Lu(end-1,:),0];
p3 = [Lu(end,:),0];
% 5.2 \textbf{while} the last 3 points in Lupper do
% not make a right turn.
c = cross(p3-p1,p2-p1);
% 6. \textbf{do} Delete the middle of the last three points from
% Lupper.
if c(3) <= 0
Lu(end-1,:) = [];
else
break
end
end
end
% 7. Put the points pn and pn?1 in a list Llower,
% with pn as the first point.
Ll = [Ps(end,:);Ps(end-1,:)];
% 8. \textbf{for} ii ¡û n?2 \textbf{downto} 1
for ii = 3:length(Ps)
% 9. \textbf{do} Append pi to Llower.
Ll = [Ll;Ps(end-ii+1,:)]; % add next point to set
% 10.1 \textbf{while} Llower contains more than 2 points
while length(Ll) > 2
p3 = [Ll(end,:),0];
p2 = [Ll(end-1,:),0];
p1 = [Ll(end-2,:),0];
% 10.2 \textbf{while} the last 3 points in Llower do
% not make a right turn.
c = cross(p3-p1,p2-p1);
% 11. \textbf{do} Delete the middle of the last 3 points from
% Llower.
if c(3) <= 0
Ll(end-1,:) = [];
else
break
end
end
end
% 12. Remove the first and the last point from Llower to
% avoid duplication of the points where the upper and lower hull meet.
Ll(1,:) = []; Ll(end,:) = [];
% 13. Append Llower to Lupper, and call the
% resulting list L .
L = [Lu;Ll];
