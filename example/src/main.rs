use bevy::{prelude::*, render::{mesh::{Indices, PrimitiveTopology, VertexAttributeValues}, render_asset::RenderAssetUsages}};
use cvoxel::CVoxels;


fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_systems(Startup, setup)
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    commands.spawn(Camera3dBundle {
        transform: Transform::from_xyz(0.0, 1.0, 2.0).looking_at(Vec3::new(0.0, 0.0, 0.0), Dir3::Y),
        ..Default::default()
    });

    commands.spawn(PointLightBundle {
        point_light: PointLight {
            intensity: 10000.0,
            color: Color::WHITE,
            ..Default::default()
        },
        transform: Transform::from_xyz(0.0, 1.0, 2.0),
        ..Default::default()
    });

    let mesh = Capsule3d::new(0.3, 0.7).mesh().build();

    let mesh_attr = mesh.attribute(Mesh::ATTRIBUTE_POSITION).unwrap();
    if let VertexAttributeValues::Float32x3(v) = mesh_attr {
        let indices = mesh.indices().unwrap();
        let dx = 0.05;
        let voxels = match indices {
            Indices::U16(ids) => CVoxels::from_indexed_mesh(&v, &ids, dx),
            Indices::U32(ids) => CVoxels::from_indexed_mesh(&v, &ids, dx)
        }.unwrap();
        let surface_vertices = voxels.surface_mesh();
        let surface_mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::default())
            .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, surface_vertices)
            .with_computed_normals();
        commands.spawn(PbrBundle {
            mesh: meshes.add(surface_mesh),
            material: materials.add(Color::WHITE),
            transform: Transform::from_xyz(0.0, 0.0, -0.0),
            ..Default::default()
        });
    }

    let mesh_handle = meshes.add(mesh);
    commands.spawn(PbrBundle {
        mesh: mesh_handle,
        material: materials.add(Color::WHITE),
        transform: Transform::from_xyz(0.0, 0.0, 0.0),
        ..Default::default()
    });
}